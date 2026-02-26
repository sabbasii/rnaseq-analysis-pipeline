#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/00_config.sh"
ensure_dirs
validate_core_inputs

need_dir "$FASTQ_DIR"

LOG="$LOG_DIR/01_qc_fastqc_multiqc.log"
echo "Logging to: $LOG"
{
  echo "REPO_ROOT=$REPO_ROOT"
  echo "FASTQ_DIR=$FASTQ_DIR"
  echo "THREADS=$THREADS"
  echo

  # FastQC on all .fastq.gz in processed_fastq
  mapfile -t FASTQS < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*.fastq.gz" | sort)
  [[ "${#FASTQS[@]}" -gt 0 ]] || { echo "ERROR: no FASTQs found in $FASTQ_DIR"; exit 1; }

  echo "FASTQ files found: ${#FASTQS[@]}"
  echo "Running FastQC..."
  fastqc \
    -o "$FASTQC_DIR" \
    -t "$THREADS" \
    --noextract \
    "${FASTQS[@]}"

  echo "Running MultiQC..."
  multiqc \
    --outdir "$MULTIQC_DIR" \
    --no-data-dir \
    --interactive \
    "$FASTQC_DIR"

  echo "Done."
} |& tee "$LOG"