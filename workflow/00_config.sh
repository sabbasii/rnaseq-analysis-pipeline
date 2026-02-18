#!/usr/bin/env bash
set -euo pipefail

# ---------------------------
# 00_config.sh
# Shared config for the pipeline.
# Source this file (do not execute).
# ---------------------------

# Resolve REPO_ROOT robustly from this file's location
# Resolve script directory robustly (bash + zsh compatible)
if [ -n "${BASH_SOURCE:-}" ]; then
  SCRIPT_PATH="${BASH_SOURCE[0]}"
elif [ -n "${(%):-%x}" ]; then
  SCRIPT_PATH="${(%):-%x}"
else
  SCRIPT_PATH="$0"
fi

SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
export REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"


# -------- Paths (repo-local) --------
export INPUT_DIR="$REPO_ROOT/input_data"
export RAW_FASTQ_DIR="$INPUT_DIR/raw_fastq"
export FASTQ_DIR="$INPUT_DIR/processed_fastq"
export MD5_DIR="$INPUT_DIR/md5"

export META_DIR="$REPO_ROOT/metadata"
export SAMPLES_TSV="$META_DIR/samples.tsv"

export REF_DIR="$REPO_ROOT/reference"
export GENOME_DIR="$REF_DIR/genome"
export STAR_INDEX_DIR="$REF_DIR/star_index"

export RESULTS_DIR="$REPO_ROOT/results"
export QC_DIR="$RESULTS_DIR/qc"
export ALIGN_DIR="$RESULTS_DIR/align"
export COUNT_DIR="$RESULTS_DIR/counts"
export LOG_DIR="$RESULTS_DIR/logs"

export FASTQC_DIR="$QC_DIR/fastqc"
export MULTIQC_DIR="$QC_DIR/multiqc"

export BAM_DIR="$ALIGN_DIR/bam"
export STAR_LOG_DIR="$ALIGN_DIR/star_logs"

# -------- Resources --------
# Override at runtime, e.g.: THREADS=32 bash workflow/02_align_star.sh
export THREADS="${THREADS:-16}"

# -------- Helpers --------
ensure_dirs() {
  mkdir -p \
    "$QC_DIR" "$ALIGN_DIR" "$COUNT_DIR" "$LOG_DIR" \
    "$FASTQC_DIR" "$MULTIQC_DIR" \
    "$BAM_DIR" "$STAR_LOG_DIR"
}

need_file() {
  local f="$1"
  [[ -f "$f" ]] || { echo "ERROR: missing file: $f" >&2; exit 1; }
}

need_dir() {
  local d="$1"
  [[ -d "$d" ]] || { echo "ERROR: missing directory: $d" >&2; exit 1; }
}

# Check core inputs exist
validate_core_inputs() {
  need_dir "$FASTQ_DIR"
  need_file "$SAMPLES_TSV"
  need_dir "$STAR_INDEX_DIR"
  # Genome + annotation (used by STAR + featureCounts)
  ls "$GENOME_DIR"/*.fna >/dev/null 2>&1 || { echo "ERROR: no .fna found in $GENOME_DIR" >&2; exit 1; }
  ls "$GENOME_DIR"/*.gtf >/dev/null 2>&1 || { echo "ERROR: no .gtf found in $GENOME_DIR" >&2; exit 1; }
}