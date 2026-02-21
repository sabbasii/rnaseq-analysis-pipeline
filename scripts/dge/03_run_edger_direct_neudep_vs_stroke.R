#!/usr/bin/env Rscript
# ------------------------------------------------------------
# scripts/dge/03_run_edger_direct_neudep_vs_stroke.R
#
# Purpose
#   Direct (cross-sectional) edgeR QL differential expression:
#     NeuDep (treated) vs Stroke at:
#       - 3h  (reperfusion)
#       - 24h (follow-up)
#
#   This is NOT a trajectory (difference-in-differences) contrast.
#   It answers: at each timepoint, is expression higher/lower in
#   NeuDep vs Stroke?
#
# Inputs (from 01_prepare_edger_inputs.R)
#   - dge/inputs/matrices/counts_matrix.tsv
#   - dge/inputs/metadata/metadata.tsv
#
# Outputs
#   - results/dge/model_direct_neudep_vs_stroke/tables/<contrast>.all.tsv
#   - results/dge/model_direct_neudep_vs_stroke/tables_filtered/<contrast>.edger_fdr25_logfc1p25.tsv
#
# Contrasts written
#   - neudep_vs_stroke_3h    = treated_3  - stroke_3
#   - neudep_vs_stroke_24h   = treated_24 - stroke_24
#
# Run (from repo root)
#   Rscript scripts/dge/03_run_edger_direct_neudep_vs_stroke.R
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

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

normalize_treatment <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  ifelse(grepl("^sham$", x0), "sham",
  ifelse(grepl("^stroke$", x0) | grepl("^mcao$", x0), "stroke",
  ifelse(grepl("anti", x0) | grepl("ly6g", x0) | grepl("neude", x0) | grepl("^treated$", x0), "treated",
         NA_character_)))
}

normalize_timepoint <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  ifelse(grepl("baseline", x0) | grepl("^bl$", x0) | grepl("^0$", x0), 0,
  ifelse(grepl("3h", x0) | grepl("reperfusion", x0) | grepl("3hpr", x0), 3,
  ifelse(grepl("24h", x0) | grepl("follow", x0) | grepl("post.?stroke", x0), 24,
         NA_integer_)))
}

repo_root <- resolve_repo_root()

counts_in <- file.path(repo_root, "dge", "inputs", "matrices", "counts_matrix.tsv")
meta_in   <- file.path(repo_root, "dge", "inputs", "metadata", "metadata.tsv")

need_file(counts_in, "counts_matrix.tsv")
need_file(meta_in,   "metadata.tsv")

model_dir  <- file.path(repo_root, "results", "dge", "model_direct_neudep_vs_stroke")
tables_dir <- file.path(model_dir, "tables")
filt_dir   <- file.path(model_dir, "tables_filtered")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filt_dir,   recursive = TRUE, showWarnings = FALSE)

counts_df <- fread(counts_in, data.table = FALSE)
meta_df   <- fread(meta_in,   data.table = FALSE)

if (!"Geneid" %in% colnames(counts_df)) stop("counts_matrix.tsv must contain Geneid column.")
if (!"SampleID" %in% colnames(meta_df)) stop("metadata.tsv must contain SampleID column.")
if (!"treatment" %in% colnames(meta_df)) stop("metadata.tsv must contain treatment column.")
if (!"timePoint" %in% colnames(meta_df)) stop("metadata.tsv must contain timePoint column.")

sample_ids <- setdiff(colnames(counts_df), "Geneid")

meta_df$SampleID <- as.character(meta_df$SampleID)
if (!all(sample_ids %in% meta_df$SampleID)) {
  missing <- setdiff(sample_ids, meta_df$SampleID)
  stop("Counts contain SampleIDs not found in metadata:\n  ", paste(missing, collapse = "\n  "))
}
meta_df <- meta_df[match(sample_ids, meta_df$SampleID), , drop = FALSE]
stopifnot(identical(meta_df$SampleID, sample_ids))

meta_df$Treatment <- normalize_treatment(meta_df$treatment)
meta_df$Time_hr   <- normalize_timepoint(meta_df$timePoint)

if (anyNA(meta_df$Treatment)) {
  bad <- unique(meta_df$treatment[is.na(meta_df$Treatment)])
  stop("Unmapped treatment values in metadata:\n  ", paste(bad, collapse = "\n  "))
}
if (anyNA(meta_df$Time_hr)) {
  bad <- unique(meta_df$timePoint[is.na(meta_df$Time_hr)])
  stop("Unmapped timePoint values in metadata:\n  ", paste(bad, collapse = "\n  "))
}

meta_df$treat_time <- factor(paste(meta_df$Treatment, meta_df$Time_hr, sep = "_"))

count_mat <- as.matrix(counts_df[, sample_ids, drop = FALSE])
rownames(count_mat) <- counts_df$Geneid
storage.mode(count_mat) <- "integer"

y <- DGEList(counts = count_mat, genes = data.frame(Geneid = counts_df$Geneid))

keep <- filterByExpr(y, group = meta_df$treat_time)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~0 + treat_time, data = meta_df)
colnames(design) <- levels(meta_df$treat_time)

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

needed_levels <- c("stroke_3", "stroke_24", "treated_3", "treated_24")
missing_lvls <- setdiff(needed_levels, colnames(design))
if (length(missing_lvls) > 0) {
  stop(
    "Design is missing required treat_time levels:\n  ",
    paste(missing_lvls, collapse = "\n  "),
    "\nNeed stroke/treated at 3h and 24h."
  )
}

contrast_mat <- makeContrasts(
  neudep_vs_stroke_3h  = treated_3  - stroke_3,
  neudep_vs_stroke_24h = treated_24 - stroke_24,
  levels = design
)

LOGFC_CUTOFF <- log2(1.25)
FDR_CUTOFF   <- 0.25

run_one <- function(name) {
  qlf <- glmQLFTest(fit, contrast = contrast_mat[, name])
  tt <- topTags(qlf, n = Inf, adjust.method = "BH")$table
  tt <- as.data.frame(tt)
  tt$Geneid <- rownames(tt)
  tt <- tt[, c("Geneid", setdiff(colnames(tt), "Geneid")), drop = FALSE]

  out_all <- file.path(tables_dir, paste0(name, ".all.tsv"))
  fwrite(tt, out_all, sep = "\t", quote = FALSE, na = "NA")

  tt_filt <- tt[
    !is.na(tt$FDR) & tt$FDR <= FDR_CUTOFF &
    !is.na(tt$logFC) & abs(tt$logFC) >= LOGFC_CUTOFF,
    , drop = FALSE
  ]
  out_filt <- file.path(filt_dir, paste0(name, ".edger_fdr25_logfc1p25.tsv"))
  fwrite(tt_filt, out_filt, sep = "\t", quote = FALSE, na = "NA")
}

invisible(lapply(colnames(contrast_mat), run_one))

cat("Done. Wrote direct contrasts to:\n  ", model_dir, "\n", sep = "")

# Run it
# Rscript scripts/dge/03_run_edger_direct_neudep_vs_stroke.R

# Output files:
# neudep_vs_stroke_3h.all.tsv
# neudep_vs_stroke_24h.all.tsv