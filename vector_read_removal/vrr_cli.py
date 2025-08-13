import argparse
import logging
import os
import sys
from typing import Optional

from .vrr_funcs import run_filter


def positive_int(value: str) -> int:
    iv = int(value)
    if iv < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return iv


def parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="vrr",
        description=(
            "Filter paired-end FASTQ reads by removing pairs that align to a reference "
            "within a given edit distance (NM) threshold using minimap2."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    req = p.add_argument_group("Required")
    req.add_argument("-1", "--r1", required=True, help="R1 FASTQ (optionally .gz)")
    req.add_argument("-2", "--r2", required=True, help="R2 FASTQ (optionally .gz)")
    req.add_argument("-r", "--ref", required=True, help="Reference FASTA")
    req.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        help=(
            "Output prefix for filtered FASTQs (default: <prefix>_R1.fastq / <prefix>_R2.fastq; "
            "use --gzip-output for .fastq.gz)"
        ),
    )

    opt = p.add_argument_group("Options")
    opt.add_argument("-t", "--threads", type=int, default=None, help="Threads for minimap2")
    opt.add_argument("--nm-threshold", type=positive_int, default=4, help="NM edit distance cutoff")
    opt.add_argument("--minimap2", default="minimap2", help="Path or name of minimap2 executable")
    opt.add_argument(
        "--extra-mm2-args",
        default=None,
        help="Extra arguments to pass to minimap2 (quoted string)",
    )
    opt.add_argument("--tmp-dir", default=None, help="Temporary directory to use")
    opt.add_argument(
        "--gzip-output",
        action="store_true",
        help="Write outputs as .fastq.gz (default writes uncompressed .fastq)",
    )
    opt.add_argument("--gzip-level", type=int, default=6, help="Compression level for outputs")

    log = p.add_argument_group("Logging")
    log.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )

    return p.parse_args(argv)


essential_files = ("r1", "r2", "ref")


def main(argv: Optional[list] = None) -> int:
    args = parse_args(argv)

    logging.basicConfig(level=getattr(logging, args.log_level), format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("vrr_cli")

    # Basic input checks
    for key in ("r1", "r2", "ref"):
        path = getattr(args, key)
        if not os.path.exists(path):
            logger.error("Input not found: %s", path)
            return 2

    try:
        out_r1, out_r2, stats = run_filter(
            r1_path=args.r1,
            r2_path=args.r2,
            ref_fasta=args.ref,
            out_prefix=args.out_prefix,
            nm_threshold=args.nm_threshold,
            threads=args.threads,
            minimap2=args.minimap2,
            extra_mm2_args=args.extra_mm2_args,
            tmp_dir=args.tmp_dir,
            gzip_output=args.gzip_output,
            gzip_level=args.gzip_level,
            logger=logging.getLogger("vector_read_removal"),
        )
    except Exception as e:
        logger.error("Failed: %s", e)
        return 1

    logger.info("Output R1: %s", out_r1)
    logger.info("Output R2: %s", out_r2)
    logger.info(
        "Stats: total_pairs=%d kept=%d dropped=%d nm_threshold=%d runtime_sec=%.3f",
        stats.get("total_pairs", 0),
        stats.get("kept_pairs", 0),
        stats.get("dropped_pairs", 0),
        stats.get("nm_threshold", 0),
        stats.get("runtime_sec", 0.0),
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
