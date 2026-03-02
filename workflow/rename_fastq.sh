#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/00_config.sh"

need_dir "$RAW_FASTQ_DIR"
mkdir -p "$FASTQ_DIR"

echo "Renaming FASTQ files..."
echo "Source: $RAW_FASTQ_DIR"
echo "Destination: $FASTQ_DIR"

for file in "$RAW_FASTQ_DIR"/*.fastq.gz; do

    base=$(basename "$file")

    # Remove extension
    name="${base%.fastq.gz}"

    # Extract BoxPosition_R1 or BoxPosition_R2
    newName="${name##*.}.fastq.gz"

    dest="$FASTQ_DIR/$newName"

    if [[ -f "$dest" ]]; then
        echo "Skipping existing file: $newName"
    else
        cp "$file" "$dest"
        echo "Created: $newName"
    fi

done

echo "Renaming complete."