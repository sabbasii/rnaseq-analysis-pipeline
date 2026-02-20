#!/usr/bin/env Rscript
# ------------------------------------------------------------
# scripts/dge/02_run_edger_main_effects.R
#
# Purpose
#   Run edgeR QL differential expression for the main analysis:
#     Treatment x TimePoint (grouped factor: treat_time)
#
# Inputs (from 01_prepare_edger_inputs.R)
#   - dge/inputs/matrices/counts_matrix.tsv
#   - dge/inputs/metadata/metadata.tsv
#
# Outputs (PRIMARY edgeR model outputs)
#   - results/dge/model_treat_time/tables/<contrast>.all.tsv
#   - results/dge/model_treat_time/tables_filtered/<contrast>.edger_fdr25_logfc1p25.tsv
#   - results/dge/model_treat_time/summary/diffGenes.edger.tsv
#
# Run (from repo root)
#   Rscript scripts/dge/02_run_edger_main_effects.R
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

model_dir   <- file.path(repo_root, "results", "dge", "model_treat_time")
tables_dir  <- file.path(model_dir, "tables")
filt_dir    <- file.path(model_dir, "tables_filtered")
sum_dir     <- file.path(model_dir, "summary")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filt_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir,    recursive = TRUE, showWarnings = FALSE)

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

needed_levels <- c("sham_0","sham_3","sham_24","stroke_0","stroke_3","stroke_24","treated_0","treated_3","treated_24")
missing_lvls <- setdiff(needed_levels, colnames(design))
if (length(missing_lvls) > 0) {
  stop(
    "Design is missing required treat_time levels:\n  ",
    paste(missing_lvls, collapse = "\n  "),
    "\nCheck that metadata includes sham/stroke/treated at time 0,3,24."
  )
}

contrast_mat <- makeContrasts(
  stroke_at_reperfusion =
    (stroke_3  - stroke_0)  - (sham_3   - sham_0),
  stroke_at_followup =
    (stroke_24 - stroke_3)  - (sham_24  - sham_3),
  treatment_at_reperfusion =
    (treated_3  - treated_0) - (stroke_3 - stroke_0),
  treatment_at_followup =
    (treated_24 - treated_3) - (stroke_24 - stroke_3),
  levels = design
)

LOGFC_CUTOFF <- log2(1.25)
FDR_CUTOFF   <- 0.25

all_tables <- list()

run_one <- function(name) {
  qlf <- glmQLFTest(fit, contrast = contrast_mat[, name])
  tt <- topTags(qlf, n = Inf, adjust.method = "BH")$table
  tt <- as.data.frame(tt)
  tt$Geneid <- rownames(tt)
  tt <- tt[, c("Geneid", setdiff(colnames(tt), "Geneid")), drop = FALSE]

  out_all <- file.path(tables_dir, paste0(name, ".all.tsv"))
  fwrite(tt, out_all, sep = "\t", quote = FALSE, na = "NA")

  if (!all(c("FDR","logFC") %in% colnames(tt)))
    stop("Expected FDR and logFC columns in edgeR output.")

  tt_filt <- tt[
    !is.na(tt$FDR) & tt$FDR <= FDR_CUTOFF &
    !is.na(tt$logFC) & abs(tt$logFC) >= LOGFC_CUTOFF,
    , drop = FALSE
  ]

  out_filt <- file.path(filt_dir, paste0(name, ".edger_fdr25_logfc1p25.tsv"))
  fwrite(tt_filt, out_filt, sep = "\t", quote = FALSE, na = "NA")

  all_tables[[name]] <<- tt
}

invisible(lapply(colnames(contrast_mat), run_one))

# Wide logFC summary (UNFILTERED)
diffGenes <- data.frame(
  Geneid = sort(unique(unlist(lapply(all_tables, function(d) d$Geneid)))),
  stringsAsFactors = FALSE
)

for (nm in names(all_tables)) {
  tmp <- all_tables[[nm]][, c("Geneid","logFC")]
  colnames(tmp) <- c("Geneid", nm)
  diffGenes <- merge(diffGenes, tmp, by = "Geneid", all.x = TRUE, sort = FALSE)
}

diffGenes <- diffGenes[
  rowSums(is.na(diffGenes[, setdiff(colnames(diffGenes), "Geneid"), drop = FALSE])) < (ncol(diffGenes) - 1),
  , drop = FALSE
]

fwrite(diffGenes, file.path(sum_dir, "diffGenes.edger.tsv"), sep = "\t", quote = FALSE, na = "NA")

cat("Done. Wrote PRIMARY edgeR outputs to:\n  ", model_dir, "\n", sep = "")
cat("Next: run interpretation/post-processing:\n  Rscript scripts/dge/07_summarize_edger_treat_time_results.R\n")