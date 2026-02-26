#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/00_config.sh"
ensure_dirs

need_dir "$STAR_LOG_DIR"

# ---- thresholds (override at runtime if you want) ----
# Example: UNIQ_MIN=60 MULTI_MAX=20 UNMAP_MAX=20 bash workflow/03_qc_postalign_star.sh
UNIQ_MIN="${UNIQ_MIN:-60}"     # flag if uniquely mapped % < this
MULTI_MAX="${MULTI_MAX:-20}"   # flag if multi-mapped % > this
UNMAP_MAX="${UNMAP_MAX:-20}"   # flag if unmapped % > this

out_dir="$RESULTS_DIR/qc/postalign"
mkdir -p "$out_dir"
out_tsv="$out_dir/star_alignment_summary.tsv"
flag_tsv="$out_dir/star_alignment_flagged.tsv"

# Header
printf "Sample\tUniquely_mapped_percent\tMulti_mapped_percent\tUnmapped_percent\tNumber_of_input_reads\tUniquely_mapped_reads\n" > "$out_tsv"

# Build summary
find "$STAR_LOG_DIR" -type f -name "*Log.final.out" | sort | while read -r f; do
  sample="$(basename "$f")"
  sample="${sample%Log.final.out}"     # keeps your existing naming prefix

  input_reads="$(awk -F '|\t' '$1 ~ /Number of input reads/ {gsub(/ /,"",$2); print $2}' "$f")"
  uniq_reads="$(awk -F '|\t' '$1 ~ /Uniquely mapped reads number/ {gsub(/ /,"",$2); print $2}' "$f")"
  uniq_pct="$(awk -F '|\t' '$1 ~ /Uniquely mapped reads %/ {gsub(/%| /,"",$2); print $2}' "$f")"
  multi_pct="$(awk -F '|\t' '$1 ~ /% of reads mapped to multiple loci/ {gsub(/%| /,"",$2); print $2}' "$f")"
  unmapped_pct="$(awk -F '|\t' '$1 ~ /% of reads unmapped/ {gsub(/%| /,"",$2); s+=$2} END{print s+0}' "$f")"

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample" "${uniq_pct:-NA}" "${multi_pct:-NA}" "${unmapped_pct:-NA}" "${input_reads:-NA}" "${uniq_reads:-NA}" \
    >> "$out_tsv"
done

# Flag + print to terminal
printf "Sample\tUniquely_mapped_percent\tMulti_mapped_percent\tUnmapped_percent\tReason\n" > "$flag_tsv"

echo
echo "Post-alignment QC thresholds:"
echo "  Uniquely mapped % < ${UNIQ_MIN}  OR  Multi-mapped % > ${MULTI_MAX}  OR  Unmapped % > ${UNMAP_MAX}"
echo

flag_count=0

# Skip header (NR==1)
tail -n +2 "$out_tsv" | while IFS=$'\t' read -r sample uniq_pct multi_pct unmapped_pct input_reads uniq_reads; do
  # Treat NA as not-checkable
  [[ "$uniq_pct" == "NA" || "$multi_pct" == "NA" || "$unmapped_pct" == "NA" ]] && continue

  reason=()

  awk_cmp_lt() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a < b)}'; }
  awk_cmp_gt() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a > b)}'; }

  if awk_cmp_lt "$uniq_pct" "$UNIQ_MIN"; then
    reason+=("uniq<${UNIQ_MIN}")
  fi
  if awk_cmp_gt "$multi_pct" "$MULTI_MAX"; then
    reason+=("multi>${MULTI_MAX}")
  fi
  if awk_cmp_gt "$unmapped_pct" "$UNMAP_MAX"; then
    reason+=("unmap>${UNMAP_MAX}")
  fi

  if [[ "${#reason[@]}" -gt 0 ]]; then
    flag_count=$((flag_count+1))
    reason_str="$(IFS=,; echo "${reason[*]}")"
    printf "%s\t%s\t%s\t%s\t%s\n" "$sample" "$uniq_pct" "$multi_pct" "$unmapped_pct" "$reason_str" >> "$flag_tsv"
    printf "FLAGGED: %s  (uniq=%s, multi=%s, unmap=%s)  reason=%s\n" \
      "$sample" "$uniq_pct" "$multi_pct" "$unmapped_pct" "$reason_str"
  fi
done

echo
echo "Wrote:"
echo "  $out_tsv"
echo "  $flag_tsv"
echo "Note: flagged samples are also printed above."


# Run
# bash workflow/03_qc_postalign.sh

# Optional: customize thresholds
# UNIQ_MIN=70 MULTI_MAX=15 UNMAP_MAX=15 bash workflow/03_qc_postalign_star.sh

# Output:
# results/qc/postalign/star_alignment_summary.tsv (all samples)
# results/qc/postalign/star_alignment_flagged.tsv (only flagged)
