#!/usr/bin/env Rscript
# ------------------------------------------------------------
# scripts/dge/07_summarize_edger_treat_time_results.R
#
# Purpose (POST-PROCESSING / INTERPRETATION)
#   Reads PRIMARY edgeR outputs produced by script 02:
#     results/dge/model_treat_time/tables/<contrast>.all.tsv
#
# It then creates TWO kinds of outputs:
#
#   A) "Significance-filtered tables" (statistics-based filtering)
#      - Uses FDR + |logFC| thresholds only
#      - Writes per-contrast filtered tables and a wide logFC table
#
#   B) "Interpretation-specific gene lists" (logic-based lists)
#      1) Delta-based rescue/depletion lists (3h and 24h):
#         delta = |treatment_logFC - stroke_logFC| >= DELTA
#      2) Additional 3h-only lists (moved from script 02):
#         - restoredFunction_3h_simple
#         - upregulated_treated_3h
#         - change_in_regulation_treated_3h
#         - downregulated_treated_3h
#         - furtherdown_treated_3h
#
#
# Outputs written under a clear subfolder:
#   results/dge/model_treat_time/postproc_delta<delta>/
#
# How to Run (from repo root)
#   Rscript scripts/dge/07_summarize_edger_treat_time_results.R

# Argument
#   --delta <numeric>  (default = 3)
#
# Example Run with delta = 2 (from repo root)
#   Rscript scripts/dge/07_summarize_edger_treat_time_results.R --delta 2
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------
# Args
# ---------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  delta <- 3
  idx <- which(args == "--delta")
  if (length(idx) > 0) {
    if (idx == length(args)) stop("--delta requires a value")
    val <- as.numeric(args[idx + 1])
    if (is.na(val) || val <= 0) stop("delta must be a positive number")
    delta <- val
  }
  list(delta = delta)
}

ARGS <- parse_args()
DELTA <- ARGS$delta
delta_tag <- paste0("delta", DELTA)

# ---------------------------
# Repo root + helpers
# ---------------------------
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

need_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
}

need_cols <- function(df, cols, label) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) stop(label, " is missing columns: ", paste(miss, collapse = ", "))
}

repo_root <- resolve_repo_root()

model_root <- file.path(repo_root, "results", "dge", "model_treat_time")
tables_dir <- file.path(model_root, "tables")

contrasts <- c(
  "stroke_at_reperfusion",
  "stroke_at_followup",
  "treatment_at_reperfusion",
  "treatment_at_followup"
)

table_paths <- file.path(tables_dir, paste0(contrasts, ".all.tsv"))
for (i in seq_along(table_paths)) need_file(table_paths[i], paste0("table for ", contrasts[i]))

# ---------------------------
# Output folders (separate!)
# ---------------------------
post_root <- file.path(model_root, paste0("postproc_", delta_tag))
filtered_dir <- file.path(post_root, "tables_filtered")
summary_dir  <- file.path(post_root, "summary")
lists_dir    <- file.path(post_root, "gene_lists")

dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(lists_dir,    recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Thresholds for "significant" tables
# ---------------------------
sig_fdr <- 0.25
fc_cut  <- log2(1.25)

read_one <- function(path) {
  df <- fread(path, data.table = FALSE)
  need_cols(df, c("Geneid", "logFC", "PValue", "FDR"), basename(path))
  df$Geneid <- as.character(df$Geneid)
  df
}

tabs <- setNames(lapply(table_paths, read_one), contrasts)

# ---------------------------
# A) Write significance-filtered tables (stats only)
# ---------------------------
write_sig <- function(name, df) {
  sig <- df[
    !is.na(df$FDR) & !is.na(df$logFC) &
      df$FDR <= sig_fdr & abs(df$logFC) >= fc_cut,
    , drop = FALSE
  ]
  out_csv <- file.path(filtered_dir, paste0(name, ".sig_fdr25_logfc1p25.", delta_tag, ".csv"))
  fwrite(sig, out_csv, sep = ",", quote = FALSE, na = "NA")
  invisible(out_csv)
}
invisible(mapply(write_sig, names(tabs), tabs))

# ---------------------------
# A) Build wide logFC table (postproc copy; unfiltered)
# ---------------------------
diffGenes <- data.frame(
  Geneid = unique(unlist(lapply(tabs, `[[`, "Geneid"))),
  stringsAsFactors = FALSE
)

add_logfc_col <- function(out_df, name, df) {
  tmp <- df[, c("Geneid", "logFC"), drop = FALSE]
  colnames(tmp)[2] <- name
  merge(out_df, tmp, by = "Geneid", all.x = TRUE, all.y = TRUE, sort = FALSE)
}

for (nm in names(tabs)) {
  diffGenes <- add_logfc_col(diffGenes, nm, tabs[[nm]])
}

diffGenes <- diffGenes[
  !apply(diffGenes[, setdiff(colnames(diffGenes), "Geneid"), drop = FALSE], 1, function(x) all(is.na(x))),
  , drop = FALSE
]

diff_out <- file.path(summary_dir, paste0("diffGenes.postproc.", delta_tag, ".csv"))
fwrite(diffGenes, diff_out, sep = ",", quote = FALSE, na = "NA")

# ---------------------------
# Helper: write gene lists
# ---------------------------
write_list <- function(vec, filename) {
  vec <- unique(na.omit(as.character(vec)))
  out <- file.path(lists_dir, paste0(filename, ".", delta_tag, ".txt"))
  fwrite(data.frame(Geneid = vec), out, sep = "\t", quote = FALSE, na = "NA", col.names = FALSE)
}

# ---------------------------
# B1) Delta-based lists (3h and 24h)
# ---------------------------
t3  <- "treatment_at_reperfusion"
s3  <- "stroke_at_reperfusion"
t24 <- "treatment_at_followup"
s24 <- "stroke_at_followup"

delta_restore <- function(treat_col, stroke_col) {
  df <- diffGenes[!is.na(diffGenes[[treat_col]]) & !is.na(diffGenes[[stroke_col]]), , drop = FALSE]
  delta <- abs(df[[treat_col]] - df[[stroke_col]])
  df$Geneid[df[[treat_col]] > 0 & df[[stroke_col]] < 0 & delta >= DELTA]
}

delta_deplete <- function(treat_col, stroke_col) {
  df <- diffGenes[!is.na(diffGenes[[treat_col]]) & !is.na(diffGenes[[stroke_col]]), , drop = FALSE]
  delta <- abs(df[[treat_col]] - df[[stroke_col]])
  df$Geneid[df[[treat_col]] < 0 & df[[stroke_col]] > 0 & delta >= DELTA]
}

write_list(delta_restore(t3, s3),  "restoredFunction_3h_delta")
write_list(delta_deplete(t3, s3),  "depletedFunction_3h_delta")
write_list(delta_restore(t24, s24), "restoredFunction_24h_delta")
write_list(delta_deplete(t24, s24), "depletedFunction_24h_delta")

# ---------------------------
# B2) Additional interpretation lists (moved from script 02)
# These are "direction-pattern" lists at 3h, without delta threshold.
# ---------------------------
if (all(c(t3, s3) %in% colnames(diffGenes))) {

  restoredFunction_3h_simple <- diffGenes$Geneid[
    !is.na(diffGenes[[t3]]) & !is.na(diffGenes[[s3]]) &
      diffGenes[[t3]] > 0 & diffGenes[[s3]] < 0
  ]

  upregulated_treated_3h <- diffGenes$Geneid[
    !is.na(diffGenes[[t3]]) &
      diffGenes[[t3]] > 0 & is.na(diffGenes[[s3]])
  ]

  change_in_regulation_treated_3h <- diffGenes$Geneid[
    !is.na(diffGenes[[t3]]) & !is.na(diffGenes[[s3]]) &
      diffGenes[[t3]] < 0 & diffGenes[[s3]] > 0
  ]

  downregulated_treated_3h <- diffGenes$Geneid[
    !is.na(diffGenes[[t3]]) &
      diffGenes[[t3]] < 0 & is.na(diffGenes[[s3]])
  ]

  furtherdown_treated_3h <- diffGenes$Geneid[
    !is.na(diffGenes[[t3]]) & !is.na(diffGenes[[s3]]) &
      (diffGenes[[t3]] < diffGenes[[s3]]) & (diffGenes[[s3]] < 0)
  ]

  write_list(restoredFunction_3h_simple,       "restoredFunction_3h_simple")
  write_list(upregulated_treated_3h,           "upregulated_treated_3h")
  write_list(change_in_regulation_treated_3h,  "change_in_regulation_treated_3h")
  write_list(downregulated_treated_3h,         "downregulated_treated_3h")
  write_list(furtherdown_treated_3h,           "furtherdown_treated_3h")
}

cat(
  "Done.\n\n",
  "PRIMARY edgeR outputs are produced by script 02 under:\n  ", model_root, "\n\n",
  "This script writes POST-PROCESSING outputs under:\n  ", post_root, "\n\n",
  "How to interpret outputs:\n",
  "  - tables_filtered/*.sig_* : stats-based filtering only (FDR + |logFC|)\n",
  "  - gene_lists/*_delta*     : interpretation lists using delta threshold (|treat - stroke| >= delta)\n",
  "  - gene_lists/*_simple*    : interpretation lists based on direction patterns (no delta)\n\n",
  "Delta threshold used: ", DELTA, "\n",
  sep = ""
)