# vector_read_removal
vector_read_removal is a pip-installable Python CLI tool to filter paired-end FASTQ reads by removing pairs that align to a reference within a given edit distance (NM) threshold using minimap2. Unaligned reads and reads with NM greater than the threshold are kept. Outputs are FASTQs (default uncompressed; pass `--gzip-output` to write `.fastq.gz`).

## Requirements
* __Python__: >= 3.9
* __Python packages__: `pysam`
* __External tool__: `minimap2` must be installed and available on your PATH

Check minimap2: `minimap2 --version`

## Installation
Install from a local clone:

```bash
pip install .
# or for development
pip install -e .
```

## Command Line Usage
The CLI is installed as `vrr`.

```bash
# default (uncompressed outputs)
vrr -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -r reference.fasta -o filtered \
  -t 8 --nm-threshold 2 --minimap2 minimap2

# write gzipped outputs
vrr -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -r reference.fasta -o filtered \
  --gzip-output --gzip-level 6
```

This produces by default:
* `filtered_R1.fastq`
* `filtered_R2.fastq`

With `--gzip-output`:
* `filtered_R1.fastq.gz`
* `filtered_R2.fastq.gz`

### Options
* `-1/--r1` R1 FASTQ (optionally `.gz`)
* `-2/--r2` R2 FASTQ (optionally `.gz`)
* `-r/--ref` Reference FASTA
* `-o/--out-prefix` Output prefix
* `-t/--threads` Threads for minimap2 (default: all available)
* `--nm-threshold` NM edit distance cutoff (default: 4)
* `--minimap2` Path or name of minimap2 executable (default: minimap2)
* `--extra-mm2-args` Extra minimap2 args (quoted string), e.g. `"-x sr"` (already defaulted)
* `--tmp-dir` Temporary directory to use
* `--gzip-output` Write outputs as `.fastq.gz` (default writes uncompressed `.fastq`)
* `--gzip-level` Gzip compression level for outputs (default: 6)
* `--log` Path to log file (default: `<out_prefix>.log`)
* `--log-level` Logging verbosity (DEBUG/INFO/WARNING/ERROR)

## Behavior
* Aligns paired FASTQs to the reference with minimap2 (`-a -x sr`).
* Computes minimal `NM` (edit distance) per read from SAM records (primary/secondary included).
* Drops the pair if either mate has minimal `NM <= threshold` (i.e., close match to reference). Otherwise, keeps the pair.
* Unmapped reads (no alignments) are kept.

## Example
```bash
vrr -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -r vectors.fasta -o sample_filtered -t 16 --nm-threshold 2
```

## Notes
* Ensure your paired FASTQs are synchronized (same order and read names).
* No samtools dependency: parsing is done directly from the SAM output via pysam.
* Large datasets: memory usage scales with number of unique read names that aligned; future versions may add streaming/ BAM modes.

## License
MIT
