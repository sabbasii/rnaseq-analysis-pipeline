#!/usr/bin/env Rscript
# ------------------------------------------------------------
# scripts/dge/01_prepare_edger_inputs.R
#
# Purpose
#   Prepare standardized edgeR inputs from:
#     1) results/counts/gene_counts.txt (featureCounts output)
#     2) metadata/samples.with_ids.tsv (metadata with SampleID)
#
# Key behaviors
#   - Reads featureCounts while preserving full BAM-path column headers
#   - Drops any COUNT columns whose BAM filename contains "_replicate"
#   - Drops any METADATA rows whose SampleID contains "_replicate"
#   - Canonicalizes the TimePoint part of SampleID so "-" vs "." mismatches won't break joins
#
# Run (from repo root)
#   Rscript scripts/dge/01_prepare_edger_inputs.R
#
# Outputs (under dge/inputs/)
#   - dge/inputs/matrices/counts_matrix.tsv
#   - dge/inputs/metadata/metadata.tsv
#   - dge/inputs/qc/library_sizes.tsv
#   - dge/inputs/qc/prep_summary.txt
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
}

resolve_repo_root <- function() {
  env_root <- Sys.getenv("REPO_ROOT")
  if (nzchar(env_root)) return(normalizePath(env_root))

  wd <- normalizePath(getwd())
  if (file.exists(file.path(wd, "README.md")) ||
      file.exists(file.path(wd, ".git")) ||
      dir.exists(file.path(wd, "results")) ||
      dir.exists(file.path(wd, "metadata"))) {
    return(wd)
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    script_dir <- dirname(normalizePath(script_path))
    return(normalizePath(file.path(script_dir, "..", "..")))
  }

  wd
}

canon_sample_id <- function(ids) {
  ids <- as.character(ids)

  split3 <- strsplit(ids, "_", fixed = TRUE)
  bad <- vapply(split3, length, integer(1)) != 3
  if (any(bad)) {
    stop(
      "Some SampleID values are not MouseID_Group_TimePoint.\nExamples:\n  ",
      paste(unique(ids[bad])[1:min(8, sum(bad))], collapse = "\n  ")
    )
  }

  mouse <- vapply(split3, `[[`, character(1), 1)
  group <- vapply(split3, `[[`, character(1), 2)
  tp    <- vapply(split3, `[[`, character(1), 3)

  tp <- gsub("[^A-Za-z0-9]+", ".", tp)
  paste(mouse, group, tp, sep = "_")
}

extract_sample_id_from_col <- function(colname) {
  base <- basename(as.character(colname))
  base <- sub("\\.bam$", "", base)
  base <- sub("(_)?Aligned\\.out$", "", base)
  base
}

find_header_line_index <- function(path) {
  x <- readLines(path, n = 400, warn = FALSE)
  idx <- which(grepl("^Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t", x))
  if (length(idx) == 0) {
    stop("Could not find the featureCounts header line starting with: Geneid<TAB>Chr<TAB>Start<TAB>End<TAB>Strand<TAB>Length")
  }
  idx[1]
}

filterNegatives <- function(dataframe) {
  returnFrame <- dataframe |>
    filter(
      rowSums(across(all_of(sampleList)) > 0) > 0
    )
  return(returnFrame)
}

repo_root <- resolve_repo_root()

featurecounts_path <- file.path(repo_root, "results", "counts", "gene_counts.txt")
metadata_path      <- file.path(repo_root, "metadata", "samples.with_ids.tsv")

assert_file_exists(featurecounts_path, "featureCounts file")
assert_file_exists(metadata_path, "metadata file")

dge_inputs_root <- file.path(repo_root, "dge", "inputs")
counts_dir      <- file.path(dge_inputs_root, "matrices")
meta_dir        <- file.path(dge_inputs_root, "metadata")
qc_dir          <- file.path(dge_inputs_root, "qc")

dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,     recursive = TRUE, showWarnings = FALSE)

counts_matrix_out <- file.path(counts_dir, "counts_matrix.tsv")
metadata_out      <- file.path(meta_dir,   "metadata.tsv")
library_sizes_out <- file.path(qc_dir,     "library_sizes.tsv")
summary_out       <- file.path(qc_dir,     "prep_summary.txt")

samples_df <- fread(metadata_path, data.table = FALSE)

required_meta_cols <- c("position", "mouseID", "treatment", "timePoint", "ageSex", "age", "sex", "sampleID")
missing_cols <- setdiff(required_meta_cols, colnames(samples_df))
if (length(missing_cols) > 0) {
  stop("metadata/samples.with_ids.tsv is missing columns: ", paste(missing_cols, collapse = ", "))
}

is_meta_replicate <- grepl("_replicate", samples_df$sampleID, fixed = TRUE)
meta_replicate_rows <- samples_df[is_meta_replicate, , drop = FALSE]
samples_df <- samples_df[!is_meta_replicate, , drop = FALSE]
n_meta_rep <- nrow(meta_replicate_rows)

samples_df$SampleID <- canon_sample_id(samples_df$sampleID)

hdr_idx <- find_header_line_index(featurecounts_path)

featurecounts_df <- fread(
  featurecounts_path,
  data.table = FALSE,
  check.names = FALSE,
  skip = hdr_idx - 1
)

featurecounts_annotation_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
if (!all(featurecounts_annotation_cols %in% colnames(featurecounts_df))) {
  stop(
    "Unexpected featureCounts format after reading.\nFound first columns: ",
    paste(colnames(featurecounts_df)[1:min(12, ncol(featurecounts_df))], collapse = ", ")
  )
}

raw_sample_columns_all <- setdiff(colnames(featurecounts_df), featurecounts_annotation_cols)
if (length(raw_sample_columns_all) == 0) stop("No sample columns found in: ", featurecounts_path)

is_replicate_col <- grepl("_replicate", basename(raw_sample_columns_all), fixed = TRUE)
replicate_columns <- raw_sample_columns_all[is_replicate_col]
raw_sample_columns <- raw_sample_columns_all[!is_replicate_col]
n_counts_rep <- length(replicate_columns)

if (length(raw_sample_columns) == 0) {
  stop("All sample columns were filtered out as replicates. Check featureCounts header / naming.")
}

sample_ids_from_counts <- vapply(raw_sample_columns, extract_sample_id_from_col, character(1))
sample_ids_from_counts <- canon_sample_id(sample_ids_from_counts)

if (anyDuplicated(sample_ids_from_counts)) {
  dup <- unique(sample_ids_from_counts[duplicated(sample_ids_from_counts)])
  stop(
    "Duplicate SampleIDs after parsing featureCounts columns:\n  ",
    paste(dup, collapse = "\n  "),
    "\n\nFix: ensure each BAM filename is unique as MouseID_Group_TimePoint (excluding replicates, which we drop)."
  )
}

counts_matrix_df <- featurecounts_df[, c("Geneid", raw_sample_columns), drop = FALSE]
colnames(counts_matrix_df) <- c("Geneid", sample_ids_from_counts)

for (sid in sample_ids_from_counts) {
  counts_matrix_df[[sid]] <- as.integer(round(as.numeric(counts_matrix_df[[sid]])))
}

if (anyDuplicated(counts_matrix_df$Geneid)) {
  counts_matrix_df <- counts_matrix_df |>
    group_by(Geneid) |>
    summarise(across(all_of(sample_ids_from_counts), ~ as.integer(sum(.x, na.rm = TRUE))), .groups = "drop") |>
    as.data.frame()
}

sampleList <- sample_ids_from_counts
geneCounts <- counts_matrix_df
geneCounts <- filterNegatives(geneCounts)
counts_matrix_df <- geneCounts

sample_ids_from_counts <- setdiff(colnames(counts_matrix_df), "Geneid")

missing_in_meta <- setdiff(sample_ids_from_counts, samples_df$SampleID)
if (length(missing_in_meta) > 0) {
  stop(
    "These SampleIDs exist in counts but not in metadata (after canonicalization):\n  ",
    paste(missing_in_meta, collapse = "\n  ")
  )
}

missing_in_counts <- setdiff(samples_df$SampleID, sample_ids_from_counts)

samples_aligned_df <- samples_df[match(sample_ids_from_counts, samples_df$SampleID), , drop = FALSE]
stopifnot(identical(samples_aligned_df$SampleID, sample_ids_from_counts))

fwrite(counts_matrix_df, counts_matrix_out, sep = "\t", quote = FALSE, na = "NA")
fwrite(samples_aligned_df, metadata_out,     sep = "\t", quote = FALSE, na = "NA")

library_sizes <- colSums(as.matrix(counts_matrix_df[, sample_ids_from_counts, drop = FALSE]))
fwrite(
  data.frame(SampleID = names(library_sizes), LibrarySize = as.numeric(library_sizes)),
  library_sizes_out,
  sep = "\t", quote = FALSE, na = "NA"
)

replicate_note <- c(
  paste0("Dropped replicate rows from metadata (n=", n_meta_rep, ")."),
  paste0("Dropped replicate columns from counts (n=", n_counts_rep, ")."),
  ""
)

if (n_counts_rep > 0) {
  replicate_note <- c(replicate_note,
                      "Dropped count columns:",
                      paste0("  ", basename(replicate_columns)),
                      "")
}

if (n_meta_rep > 0) {
  replicate_note <- c(replicate_note,
                      "Dropped metadata rows (sampleID):",
                      paste0("  ", meta_replicate_rows$sampleID),
                      "")
}

summary_lines <- c(
  paste0("REPO_ROOT: ", repo_root),
  paste0("FeatureCounts input: ", featurecounts_path),
  paste0("Metadata input: ", metadata_path),
  paste0("Header line detected at: ", hdr_idx),
  "",
  replicate_note,
  paste0("Genes (rows): ", nrow(counts_matrix_df)),
  paste0("Samples (cols): ", length(sample_ids_from_counts)),
  "",
  paste0("Wrote counts matrix: ", counts_matrix_out),
  paste0("Wrote aligned metadata: ", metadata_out),
  paste0("Wrote library sizes: ", library_sizes_out),
  paste0("Wrote summary: ", summary_out),
  ""
)

if (length(missing_in_counts) > 0) {
  summary_lines <- c(
    summary_lines,
    "NOTE: metadata SampleIDs not found in counts (ignored):",
    paste(missing_in_counts, collapse = ", "),
    ""
  )
}

writeLines(summary_lines, summary_out)
cat(paste(summary_lines, collapse = "\n"), "\n")
