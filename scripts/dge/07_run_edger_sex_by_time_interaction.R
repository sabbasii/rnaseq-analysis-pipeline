#!/usr/bin/env Rscript
# ------------------------------------------------------------
# scripts/dge/07_run_edger_sex_by_time_interaction.R
#
# Purpose
#   Run edgeR QL differential expression for SEX effects by time:
#     sex_time = Sex_Time
#
# Inputs (from 01_prepare_edger_inputs.R)
#   - dge/inputs/matrices/counts_matrix.tsv
#   - dge/inputs/metadata/metadata.tsv
#
# Output (minimal)
#   - results/dge/model_sex_time/tables/<contrast>.all.tsv
#   - results/dge/model_sex_time/gene_lists/<list>.tsv   (Top genes, PValue <= 0.05)
#
# Run (from repo root)
#   Rscript scripts/dge/07_run_edger_sex_by_time_interaction.R
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

# ---------- Helpers ----------
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

pick_col <- function(df, candidates, label) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop("metadata.tsv is missing required column for ", label, ". Tried: ",
         paste(candidates, collapse = ", "))
  }
  hit[1]
}

normalize_timepoint <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  ifelse(grepl("baseline", x0) | grepl("^bl$", x0) | grepl("^0$", x0), 0,
  ifelse(grepl("3h", x0) | grepl("reperfusion", x0) | grepl("3hpr", x0), 3,
  ifelse(grepl("24h", x0) | grepl("follow", x0) | grepl("post.?stroke", x0), 24,
         NA_integer_)))
}

normalize_sex <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  ifelse(x0 %in% c("m","male"), "male",
  ifelse(x0 %in% c("f","female"), "female",
         NA_character_))
}

# ---------- Paths ----------
repo_root <- resolve_repo_root()

counts_in <- file.path(repo_root, "dge", "inputs", "matrices", "counts_matrix.tsv")
meta_in   <- file.path(repo_root, "dge", "inputs", "metadata", "metadata.tsv")

need_file(counts_in, "counts_matrix.tsv")
need_file(meta_in,   "metadata.tsv")

tables_dir     <- file.path(repo_root, "results", "dge", "model_sex_time", "tables")
gene_lists_dir <- file.path(repo_root, "results", "dge", "model_sex_time", "gene_lists")
dir.create(tables_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(gene_lists_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Load inputs ----------
counts_df <- fread(counts_in, data.table = FALSE)
meta_df   <- fread(meta_in,   data.table = FALSE)

if (!"Geneid" %in% colnames(counts_df)) stop("counts_matrix.tsv must contain Geneid column.")
if (!"SampleID" %in% colnames(meta_df)) stop("metadata.tsv must contain SampleID column.")

sample_ids <- setdiff(colnames(counts_df), "Geneid")

# Align metadata rows to counts column order
meta_df$SampleID <- as.character(meta_df$SampleID)
if (!all(sample_ids %in% meta_df$SampleID)) {
  missing <- setdiff(sample_ids, meta_df$SampleID)
  stop("Counts contain SampleIDs not found in metadata:\n  ", paste(missing, collapse = "\n  "))
}
meta_df <- meta_df[match(sample_ids, meta_df$SampleID), , drop = FALSE]
stopifnot(identical(meta_df$SampleID, sample_ids))

sex_col <- pick_col(meta_df, c("Sex", "sex"), "Sex")
tp_col  <- pick_col(meta_df, c("TimePoint", "timePoint"), "TimePoint")

meta_df$Sex2     <- normalize_sex(meta_df[[sex_col]])
meta_df$Time_hr  <- normalize_timepoint(meta_df[[tp_col]])

if (anyNA(meta_df$Sex2)) {
  bad <- unique(meta_df[[sex_col]][is.na(meta_df$Sex2)])
  stop("Unmapped Sex values:\n  ", paste(bad, collapse = "\n  "))
}
if (anyNA(meta_df$Time_hr)) {
  bad <- unique(meta_df[[tp_col]][is.na(meta_df$Time_hr)])
  stop("Unmapped TimePoint values:\n  ", paste(bad, collapse = "\n  "))
}

meta_df$sex_time <- factor(paste(meta_df$Sex2, meta_df$Time_hr, sep = "_"))

# ---------- Build DGEList ----------
count_mat <- as.matrix(counts_df[, sample_ids, drop = FALSE])
rownames(count_mat) <- counts_df$Geneid
storage.mode(count_mat) <- "integer"

y <- DGEList(counts = count_mat, genes = data.frame(Geneid = counts_df$Geneid))

keep <- filterByExpr(y, group = meta_df$sex_time)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~0 + sex_time, data = meta_df)
colnames(design) <- levels(meta_df$sex_time)

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# ---------- Contrasts ----------
needed_levels <- c("male_3","female_3","male_24","female_24")
missing_lvls <- setdiff(needed_levels, colnames(design))
if (length(missing_lvls) > 0) {
  stop(
    "Design is missing required sex_time levels:\n  ",
    paste(missing_lvls, collapse = "\n  "),
    "\nCheck that metadata includes both sexes at 3h and 24h."
  )
}

contrast_mat <- makeContrasts(
  sex_at_reperfusion = (male_3  - female_3),
  sex_at_followup    = (male_24 - female_24),
  levels = design
)

# ---------- Run + write one file per contrast ----------
run_one <- function(name) {
  qlf <- glmQLFTest(fit, contrast = contrast_mat[, name])
  tt <- topTags(qlf, n = Inf, adjust.method = "BH")$table
  tt <- as.data.frame(tt)
  tt$Geneid <- rownames(tt)
  tt <- tt[, c("Geneid", setdiff(colnames(tt), "Geneid")), drop = FALSE]

  out_path <- file.path(tables_dir, paste0(name, ".all.tsv"))
  fwrite(tt, out_path, sep = "\t", quote = FALSE, na = "NA")
  tt
}

res_list <- lapply(colnames(contrast_mat), run_one)
names(res_list) <- colnames(contrast_mat)

# ---------- Top genes lists (PValue <= 0.05; sorted by |logFC|) ----------
write_top_list <- function(df, out_name) {
  if (!("PValue" %in% colnames(df))) stop("Expected PValue column in edgeR table.")
  if (!("logFC"  %in% colnames(df))) stop("Expected logFC column in edgeR table.")

  out <- df[df$PValue <= 0.05, , drop = FALSE]
  if (nrow(out) > 0) out <- out[order(-abs(out$logFC)), , drop = FALSE]

  out_path <- file.path(gene_lists_dir, out_name)
  fwrite(out, out_path, sep = "\t", quote = FALSE, na = "NA")
}

write_top_list(res_list[["sex_at_reperfusion"]], "Male-Female_at_Reperfusion.tsv")
write_top_list(res_list[["sex_at_followup"]],    "Male-Female_at_Followup.tsv")

cat("Done.\nWrote tables to:\n  ", tables_dir,
    "\nWrote gene lists to:\n  ", gene_lists_dir, "\n", sep = "")
