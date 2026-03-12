#!/usr/bin/env Rscript
# ============================================================
# scripts/dge/pde/02_pde_neudep_vs_stroke_direct_summary.R
#
# Purpose
#   Summarize PDE family differential expression using DIRECT
#   NeuDep vs Stroke contrasts at:
#
#       - 3h   (neudep_vs_stroke_3h)
#       - 24h  (neudep_vs_stroke_24h)
#
#   This script reads existing edgeR results and extracts only
#   PDE genes, marking significance using either:
#
#       --fdr  <threshold>
#   OR  --pval <threshold>
#
#   (mutually exclusive; default is --fdr 0.10)
#
# Inputs
#   results/dge/model_direct_neudep_vs_stroke/tables/
#       neudep_vs_stroke_3h.all.tsv
#       neudep_vs_stroke_24h.all.tsv
#
#   metadata/gene_sets/pde_pathway_genes.txt
#   or
#   metadata/gene_sets/pde_genes.txt
#
# Output
#   results/dge/pde/neudep_vs_stroke_direct_pde_summary.tsv
#
# Significance definition
#   statistic threshold (FDR or PValue) AND
#   |logFC| >= log2(fc_cut)
#
# Defaults
#   --fdr 0.10
#   --fc  1.20
#
# Example
#   Rscript scripts/dge/pde/02_pde_neudep_vs_stroke_direct_summary.R
#
#   Rscript scripts/dge/pde/02_pde_neudep_vs_stroke_direct_summary.R \
#       --pval 0.05 --fc 1.25
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) stop(flag, " requires a value")
  args[i + 1]
}

fdr_cut  <- get_arg("--fdr", NA)
pval_cut <- get_arg("--pval", NA)
fc_cut   <- as.numeric(get_arg("--fc", "1.20"))

if (!is.na(fdr_cut))  fdr_cut  <- as.numeric(fdr_cut)
if (!is.na(pval_cut)) pval_cut <- as.numeric(pval_cut)

if (!is.na(fdr_cut) && !is.na(pval_cut)) {
  stop("Use only one of --fdr or --pval (not both).")
}

if (is.na(fdr_cut) && is.na(pval_cut)) {
  fdr_cut <- 0.10
}

if (is.na(fc_cut) || fc_cut <= 1) {
  stop("--fc must be > 1")
}

logfc_thr <- log2(fc_cut)

use_stat <- if (!is.na(pval_cut)) "PValue" else "FDR"
stat_cut <- if (!is.na(pval_cut)) pval_cut else fdr_cut

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

repo_root <- Sys.getenv("REPO_ROOT")
if (repo_root == "") repo_root <- normalizePath(getwd())

tables_dir <- file.path(
  repo_root,
  "results", "dge", "model_direct_neudep_vs_stroke", "tables"
)

tab_3h  <- file.path(tables_dir, "neudep_vs_stroke_3h.all.tsv")
tab_24h <- file.path(tables_dir, "neudep_vs_stroke_24h.all.tsv")

pde_file <- file.path(
  repo_root,
  "metadata", "gene_sets", "pde_genes.txt"
)

out_dir <- file.path(repo_root, "results", "dge", "pde")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(
  out_dir,
  "neudep_vs_stroke_direct_pde_summary.tsv"
)

need_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

need_file(tab_3h)
need_file(tab_24h)
need_file(pde_file)

# ------------------------------------------------------------
# Load PDE genes
# ------------------------------------------------------------

pde <- fread(pde_file, header = FALSE, data.table = FALSE)[[1]]
pde <- trimws(as.character(pde))
pde <- pde[nzchar(pde)]
pde <- pde[!grepl("^#", pde)]
pde <- unique(pde)

# ------------------------------------------------------------
# Load edgeR tables
# ------------------------------------------------------------

read_edger <- function(path, suffix) {

  df <- fread(path, data.table = FALSE)

  required <- c("Geneid", "logFC", "PValue", "FDR")
  miss <- setdiff(required, colnames(df))

  if (length(miss) > 0)
    stop("Missing columns in ", path, ": ", paste(miss, collapse = ", "))

  df <- df[, required]

  colnames(df) <- c(
    "Geneid",
    paste0("logFC_", suffix),
    paste0("PValue_", suffix),
    paste0("FDR_", suffix)
  )

  df
}

df3  <- read_edger(tab_3h,  "3h")
df24 <- read_edger(tab_24h, "24h")

# ------------------------------------------------------------
# Merge and filter PDE genes
# ------------------------------------------------------------

out <- merge(df3, df24, by = "Geneid", all = TRUE)

out <- out[out$Geneid %in% pde, ]

# ------------------------------------------------------------
# Significance flagging
# ------------------------------------------------------------

flag_sig <- function(logFC, PValue, FDR) {

  stat_ok <- if (use_stat == "PValue") {
    !is.na(PValue) & PValue <= pval_cut
  } else {
    !is.na(FDR) & FDR <= fdr_cut
  }

  fc_ok <- !is.na(logFC) & abs(logFC) >= logfc_thr

  stat_ok & fc_ok
}

out$Sig_3h  <- flag_sig(out$logFC_3h,  out$PValue_3h,  out$FDR_3h)
out$Sig_24h <- flag_sig(out$logFC_24h, out$PValue_24h, out$FDR_24h)

direction <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x > 0, "Up",
    ifelse(x < 0, "Down", "0"))
  )
}

out$Dir_3h  <- direction(out$logFC_3h)
out$Dir_24h <- direction(out$logFC_24h)

# ------------------------------------------------------------
# Sort (significant first)
# ------------------------------------------------------------

out$AnySig <- out$Sig_3h | out$Sig_24h

out <- out[
  order(
    -out$AnySig,
    -abs(out$logFC_3h),
    -abs(out$logFC_24h),
    out$Geneid
  ),
]

out$AnySig <- NULL

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fwrite(
  out,
  out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "Saved:\n  ", out_file, "\n",
  "Statistic used: ", use_stat, "\n",
  "Threshold: ", stat_cut, "\n",
  "|logFC| cutoff: ", round(logfc_thr, 3), "\n",
  sep = ""
)

# Run it:
# Rscript scripts/dge/pde/02_pde_neudep_vs_stroke_direct_summary.R

# Use P-value instead of FDR:
# Rscript scripts/dge/pde/02_pde_neudep_vs_stroke_direct_summary.R \
#   --pval 0.05 --fc 1.20

# Stricter FDR:
# Rscript scripts/dge/pde/02_pde_neudep_vs_stroke_direct_summary.R \
#   --fdr 0.05 --fc 1.25