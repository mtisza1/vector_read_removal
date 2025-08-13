import gzip
import contextlib
import logging
import os
import shlex
import shutil
import subprocess
import tempfile
import time
from typing import Dict, Optional, Tuple

import pysam


def _check_minimap2(minimap2: str) -> str:
    """Return the resolved path to minimap2 or raise if not found."""
    path = shutil.which(minimap2)
    if path is None:
        raise FileNotFoundError(
            f"Could not find minimap2 executable '{minimap2}'. Ensure minimap2 is installed and on PATH."
        )
    return path


def _run_minimap2(
    r1_path: str,
    r2_path: str,
    ref_fasta: str,
    threads: int,
    minimap2: str,
    extra_mm2_args: Optional[str],
    sam_out_path: str,
    logger: logging.Logger,
) -> None:
    mm2_path = _check_minimap2(minimap2)

    cmd = [mm2_path, "-a", "-x", "sr", "-t", str(threads)]
    if extra_mm2_args:
        cmd.extend(shlex.split(extra_mm2_args))
    cmd.extend([ref_fasta, r1_path, r2_path])

    logger.info("Running minimap2: %s", " ".join(shlex.quote(c) for c in cmd))

    with open(sam_out_path, "wb") as sam_fh:
        proc = subprocess.run(
            cmd,
            stdout=sam_fh,
            stderr=subprocess.PIPE,
            check=False,
        )
    if proc.returncode != 0:
        stderr = proc.stderr.decode("utf-8", errors="replace")
        raise RuntimeError(f"minimap2 failed with code {proc.returncode}:\n{stderr}")


def run_filter(
    r1_path: str,
    r2_path: str,
    ref_fasta: str,
    out_prefix: str,
    nm_threshold: int = 2,
    threads: Optional[int] = None,
    minimap2: str = "minimap2",
    extra_mm2_args: Optional[str] = None,
    tmp_dir: Optional[str] = None,
    keep_either: bool = False,
    gzip_output: bool = False,
    gzip_level: int = 6,
    logger: Optional[logging.Logger] = None,
) -> Tuple[str, str, Dict[str, int]]:
    """
    Align paired FASTQ reads against a reference with minimap2 and keep pairs that do not
    align within the NM edit distance threshold.

    Returns (out_r1_path, out_r2_path, stats_dict)
    """
    start_t = time.time()

    # Logger setup
    logger = logger or logging.getLogger("vector_read_removal")
    if not logger.handlers:
        logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    # Validate inputs
    for p in (r1_path, r2_path, ref_fasta):
        if not os.path.exists(p):
            raise FileNotFoundError(f"Input file not found: {p}")

    threads = threads or (os.cpu_count() or 1)

    # Prepare temp SAM path
    cleanup_tmpdir = False
    if tmp_dir is None:
        tmp_dir = tempfile.mkdtemp(prefix="vrr_")
        cleanup_tmpdir = True
    os.makedirs(tmp_dir, exist_ok=True)
    sam_path = os.path.join(tmp_dir, "minimap2_out.sam")

    try:
        # Run minimap2 to produce SAM
        _run_minimap2(
            r1_path=r1_path,
            r2_path=r2_path,
            ref_fasta=ref_fasta,
            threads=threads,
            minimap2=minimap2,
            extra_mm2_args=extra_mm2_args,
            sam_out_path=sam_path,
            logger=logger,
        )

        # Parse SAM and collect minimal NM per read for R1 and R2
        min_nm_r1: Dict[str, int] = {}
        min_nm_r2: Dict[str, int] = {}
        total_sam_records = 0
        with pysam.AlignmentFile(sam_path, "r") as sam:
            for aln in sam:
                total_sam_records += 1
                if aln.is_unmapped:
                    continue
                try:
                    nm = aln.get_tag("NM")
                except KeyError:
                    # No NM tag present; skip from drop consideration
                    continue
                name = aln.query_name
                if aln.is_read1:
                    prev = min_nm_r1.get(name)
                    if prev is None or nm < prev:
                        min_nm_r1[name] = nm
                elif aln.is_read2:
                    prev = min_nm_r2.get(name)
                    if prev is None or nm < prev:
                        min_nm_r2[name] = nm
                else:
                    # In rare cases, orientation flags might not indicate read1/read2; ignore
                    continue

        drop_r1 = {n for n, nm in min_nm_r1.items() if nm <= nm_threshold}
        drop_r2 = {n for n, nm in min_nm_r2.items() if nm <= nm_threshold}

        logger.info(
            "Collected NM: R1 candidates=%d, R2 candidates=%d; drop<=%d -> R1=%d, R2=%d",
            len(min_nm_r1), len(min_nm_r2), nm_threshold, len(drop_r1), len(drop_r2)
        )

        # Prepare outputs (default: uncompressed .fastq; if gzip_output -> .fastq.gz)
        if gzip_output:
            out_r1 = f"{out_prefix}_R1.fastq.gz"
            out_r2 = f"{out_prefix}_R2.fastq.gz"
        else:
            out_r1 = f"{out_prefix}_R1.fastq"
            out_r2 = f"{out_prefix}_R2.fastq"

        total_pairs = 0
        kept_pairs = 0
        dropped_pairs = 0

        with contextlib.ExitStack() as stack:
            fh1 = stack.enter_context(pysam.FastxFile(r1_path))
            fh2 = stack.enter_context(pysam.FastxFile(r2_path))
            if gzip_output:
                w1 = stack.enter_context(gzip.open(out_r1, "wt", compresslevel=gzip_level))
                w2 = stack.enter_context(gzip.open(out_r2, "wt", compresslevel=gzip_level))
            else:
                w1 = stack.enter_context(open(out_r1, "w"))
                w2 = stack.enter_context(open(out_r2, "w"))

            for rec1, rec2 in zip(fh1, fh2):
                total_pairs += 1
                name1 = rec1.name
                name2 = rec2.name

                # Default: keep only if both mates are NOT in drop sets
                if keep_either:
                    keep = (name1 not in drop_r1) or (name2 not in drop_r2)
                else:
                    keep = (name1 not in drop_r1) and (name2 not in drop_r2)

                if keep:
                    kept_pairs += 1
                    # Write R1
                    if rec1.comment:
                        header1 = f"@{rec1.name} {rec1.comment}\n"
                    else:
                        header1 = f"@{rec1.name}\n"
                    w1.write(header1)
                    w1.write(f"{rec1.sequence}\n+\n{rec1.quality}\n")

                    # Write R2
                    if rec2.comment:
                        header2 = f"@{rec2.name} {rec2.comment}\n"
                    else:
                        header2 = f"@{rec2.name}\n"
                    w2.write(header2)
                    w2.write(f"{rec2.sequence}\n+\n{rec2.quality}\n")
                else:
                    dropped_pairs += 1

        stats = {
            "total_pairs": total_pairs,
            "kept_pairs": kept_pairs,
            "dropped_pairs": dropped_pairs,
            "nm_threshold": nm_threshold,
            "threads": threads,
            "total_sam_records": total_sam_records,
            "runtime_sec": round(time.time() - start_t, 3),
        }

        logger.info(
            "Done. total_pairs=%d kept=%d dropped=%d runtime=%.3fs",
            total_pairs, kept_pairs, dropped_pairs, stats["runtime_sec"]
        )

        return out_r1, out_r2, stats

    finally:
        # Cleanup temp dir if we created it
        try:
            if cleanup_tmpdir and os.path.isdir(tmp_dir):
                # Remove created files
                for fname in os.listdir(tmp_dir):
                    try:
                        os.remove(os.path.join(tmp_dir, fname))
                    except OSError:
                        pass
                try:
                    os.rmdir(tmp_dir)
                except OSError:
                    pass
        except Exception:
            # Best-effort cleanup only
            pass
