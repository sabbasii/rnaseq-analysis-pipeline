#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/00_config.sh"
ensure_dirs
validate_core_inputs

need_file "$SAMPLES_TSV"
need_dir "$FASTQ_DIR"
need_dir "$STAR_INDEX_DIR"

LOG="$LOG_DIR/02_align_star.log"
echo "Logging to: $LOG"

{
  echo "REPO_ROOT=$REPO_ROOT"
  echo "SAMPLES_TSV=$SAMPLES_TSV"
  echo "FASTQ_DIR=$FASTQ_DIR"
  echo "STAR_INDEX_DIR=$STAR_INDEX_DIR"
  echo "THREADS=$THREADS"
  echo

  missing=0
  total=0

  # Read TSV: BoxPosition  MouseID  Group  TimePoint
  while IFS=$'\t' read -r box mouse group timepoint; do
    # skip empty lines
    [[ -n "${box:-}" ]] || continue
    total=$((total+1))

    r1="$FASTQ_DIR/${box}_R1.fastq.gz"
    r2="$FASTQ_DIR/${box}_R2.fastq.gz"

    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
      echo "MISSING FASTQ PAIR for $box  ($r1 / $r2)"
      missing=$((missing+1))
      continue
    fi

    # Sanitize timepoint/group for filenames
    safe_group="${group// /_}"
    safe_tp="${timepoint// /_}"

    prefix="$STAR_LOG_DIR/${mouse}_${safe_group}_${safe_tp}_"

    echo "STAR: $box  -> ${mouse}_${safe_group}_${safe_tp}"

    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$STAR_INDEX_DIR" \
      --readFilesIn "$r1" "$r2" \
      --readFilesCommand zcat \
      --outFileNamePrefix "$prefix" \
      --outSAMtype BAM Unsorted

    # Move BAM and rename according to documented convention
    bam_src="${prefix}Aligned.out.bam"
    bam_dst="$BAM_DIR/${mouse}_${safe_group}_${safe_tp}_Aligned.out.bam"

    if [[ -f "$bam_src" ]]; then
      mv "$bam_src" "$bam_dst"
      echo "  BAM -> $bam_dst"
    else
      echo "ERROR: STAR did not produce expected BAM: $bam_src"
      exit 1
    fi

  done < "$SAMPLES_TSV"

  echo
  echo "Samples in metadata: $total"
  echo "Missing FASTQ pairs: $missing"

  if [[ "$missing" -gt 0 ]]; then
    echo "ERROR: Some FASTQ pairs were missing. Fix processed_fastq/ then re-run."
    exit 1
  fi

  echo "Done."
} |& tee "$LOG"