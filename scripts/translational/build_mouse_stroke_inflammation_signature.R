#!/usr/bin/env Rscript
# ============================================================
# scripts/translational/build_mouse_stroke_inflammation_signature.R
#
# Purpose
#   Build a mouse-derived, neutrophil-relevant “stroke inflammation” signature
#   using effect-pattern logic across *direct* contrasts.
#
#   Step 1 (Stroke-induced genes):
#     Select genes significantly changed in Stroke vs Sham at 3h and/or 24h
#     using default thresholds:
#       - PValue <= 0.05
#       - |logFC| >= 1.2   (edgeR logFC is log2 fold-change)
#
#   Step 2 (Neutrophil-dependent prioritization):
#     Identify genes significantly changed in NeuDep vs Stroke at matched timepoints,
#     using the same default thresholds unless overridden.
#
# Inputs
#   Direct edgeR tables (all genes):
#     - results/dge/model_direct_stroke_vs_sham/tables/stroke_vs_sham_3h.all.tsv
#     - results/dge/model_direct_stroke_vs_sham/tables/stroke_vs_sham_24h.all.tsv
#     - results/dge/model_direct_neudep_vs_stroke/tables/neudep_vs_stroke_3h.all.tsv
#     - results/dge/model_direct_neudep_vs_stroke/tables/neudep_vs_stroke_24h.all.tsv
#
# Outputs
#   - metadata/gene_sets/stroke_inflammation_signature_mouse.txt
#   - results/dge/target_gene_sets/signatures/stroke_inflammation_signature.tsv
#
# NOTE ON OUTPUTS
#   This script produces two complementary outputs:
#
#   1) metadata/gene_sets/stroke_inflammation_signature_mouse.txt
#      → A clean, reusable gene list for downstream analyses
#        (enrichment, highlighting, cross-species mapping, etc.).
#
#   2) results/dge/target_gene_sets/signatures/stroke_inflammation_signature.tsv
#      → A full evidence table containing statistics supporting inclusion
#        of each gene in the signature (for transparency and reproducibility).
#
# Run examples
#   Rscript scripts/translational/build_mouse_stroke_inflammation_signature.R
#
#   Rscript scripts/translational/build_mouse_stroke_inflammation_signature.R \
#     --use_stat FDR
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------
# Args
# -----------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) stop(flag, " requires a value")
  args[i + 1]
}

# Choose which significance statistic to use (default: PValue)
use_stat <- toupper(get_arg("--use_stat", "PVALUE"))  # PVALUE | FDR
if (!use_stat %in% c("PVALUE", "FDR")) stop("--use_stat must be PVALUE or FDR.")

# New preferred flags (edgeR logFC is log2 fold-change)
p_stroke     <- as.numeric(get_arg("--p_stroke",     "0.05"))
logfc_stroke <- as.numeric(get_arg("--logfc_stroke", "1.2"))

p_neudep     <- as.numeric(get_arg("--p_neudep",     "0.05"))
logfc_neudep <- as.numeric(get_arg("--logfc_neudep", "1.2"))

# Backward-compatible flags (if provided, override)
fdr_stroke <- as.numeric(get_arg("--fdr_stroke", NA))
fc_stroke  <- as.numeric(get_arg("--fc_stroke",  NA))
fdr_neudep <- as.numeric(get_arg("--fdr_neudep", NA))
fc_neudep  <- as.numeric(get_arg("--fc_neudep",  NA))

if (!is.na(fdr_stroke)) { use_stat <- "FDR"; p_stroke <- fdr_stroke }
if (!is.na(fdr_neudep)) { use_stat <- "FDR"; p_neudep <- fdr_neudep }

if (!is.na(fc_stroke))  logfc_stroke <- log2(fc_stroke)
if (!is.na(fc_neudep))  logfc_neudep <- log2(fc_neudep)

# Timepoint inclusion (default: both)
use_3h  <- tolower(get_arg("--use_3h",  "true")) %in% c("true","t","1","yes","y")
use_24h <- tolower(get_arg("--use_24h", "true")) %in% c("true","t","1","yes","y")

# Validate
if (is.na(p_stroke)  || p_stroke  <= 0 || p_stroke  > 1) stop("Stroke p/FDR threshold must be in (0,1].")
if (is.na(p_neudep)  || p_neudep  <= 0 || p_neudep  > 1) stop("NeuDep p/FDR threshold must be in (0,1].")
if (is.na(logfc_stroke) || logfc_stroke <= 0) stop("--logfc_stroke must be > 0.")
if (is.na(logfc_neudep) || logfc_neudep <= 0) stop("--logfc_neudep must be > 0.")
if (!use_3h && !use_24h) stop("At least one of --use_3h or --use_24h must be true.")

# -----------------------
# Repo root + paths
# -----------------------
repo_root <- Sys.getenv("REPO_ROOT")
if (!nzchar(repo_root)) repo_root <- normalizePath(getwd())

need_file <- function(path, label = path) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
}

stroke_dir <- file.path(repo_root, "results", "dge", "model_direct_stroke_vs_sham", "tables")
neudep_dir <- file.path(repo_root, "results", "dge", "model_direct_neudep_vs_stroke", "tables")

stroke_3h_file  <- file.path(stroke_dir, "stroke_vs_sham_3h.all.tsv")
stroke_24h_file <- file.path(stroke_dir, "stroke_vs_sham_24h.all.tsv")
neudep_3h_file  <- file.path(neudep_dir, "neudep_vs_stroke_3h.all.tsv")
neudep_24h_file <- file.path(neudep_dir, "neudep_vs_stroke_24h.all.tsv")

if (use_3h) {
  need_file(stroke_3h_file,  "stroke_vs_sham_3h.all.tsv")
  need_file(neudep_3h_file,  "neudep_vs_stroke_3h.all.tsv")
}
if (use_24h) {
  need_file(stroke_24h_file, "stroke_vs_sham_24h.all.tsv")
  need_file(neudep_24h_file, "neudep_vs_stroke_24h.all.tsv")
}

# Outputs
sig_dir <- file.path(repo_root, "results", "dge", "target_gene_sets", "signatures")
dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)

mouse_out <- file.path(repo_root, "metadata", "gene_sets", "stroke_inflammation_signature_mouse.txt")
dir.create(dirname(mouse_out), recursive = TRUE, showWarnings = FALSE)

tsv_out <- file.path(sig_dir, "stroke_inflammation_signature.tsv")

# -----------------------
# Load + standardize tables
# -----------------------
read_edger <- function(path) {
  df <- fread(path, data.table = FALSE)
  req <- c("Geneid","logFC","PValue","FDR")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Table missing columns (", basename(path), "): ", paste(miss, collapse = ", "))
  df$Geneid <- as.character(df$Geneid)
  df
}

keep_sig <- function(stat, logfc, stat_thr, logfc_thr) {
  !is.na(stat) & stat <= stat_thr & !is.na(logfc) & abs(logfc) >= logfc_thr
}

make_pair <- function(stroke_df, neudep_df, time_label) {
  s <- stroke_df[, c("Geneid","logFC","PValue","FDR")]
  n <- neudep_df[, c("Geneid","logFC","PValue","FDR")]
  colnames(s) <- c("Geneid","logFC_stroke","P_stroke","FDR_stroke")
  colnames(n) <- c("Geneid","logFC_neudep","P_neudep","FDR_neudep")

  # Inner join: only genes present in BOTH tables are evaluated for paired inclusion
  m <- merge(s, n, by = "Geneid", all = FALSE, sort = FALSE)
  m$Time <- time_label

  stat_stroke <- if (use_stat == "FDR") m$FDR_stroke else m$P_stroke
  stat_neudep <- if (use_stat == "FDR") m$FDR_neudep else m$P_neudep

  # Per-contrast significance + effect filters
  m$Stroke_sig <- keep_sig(stat_stroke, m$logFC_stroke, p_stroke, logfc_stroke)
  m$NeuDep_sig <- keep_sig(stat_neudep, m$logFC_neudep, p_neudep, logfc_neudep)

  # Keep = union of Stroke_sig and NeuDep_sig
  m$Keep <- m$Stroke_sig | m$NeuDep_sig

  # One label per row (no duplicates)
  m$Class <- ifelse(m$Stroke_sig, "Stroke_sig",
                   ifelse(m$NeuDep_sig, "NeuDep_sig", NA))

  m
}

# -----------------------
# Build pairs
# -----------------------
pairs <- list()

if (use_3h) {
  s3 <- read_edger(stroke_3h_file)
  n3 <- read_edger(neudep_3h_file)
  pairs[["3h"]] <- make_pair(s3, n3, "3h")
}
if (use_24h) {
  s24 <- read_edger(stroke_24h_file)
  n24 <- read_edger(neudep_24h_file)
  pairs[["24h"]] <- make_pair(s24, n24, "24h")
}

all_pairs <- rbindlist(pairs, use.names = TRUE, fill = TRUE)

# -----------------------
# Build signature lists
# -----------------------
sig_table <- all_pairs[all_pairs$Keep, , drop = FALSE]
mouse_sig <- sort(unique(sig_table$Geneid))

writeLines(mouse_sig, mouse_out)
fwrite(sig_table, tsv_out, sep = "\t", quote = FALSE, na = "NA")

# -----------------------
# Diagnostic counts
# -----------------------
cat("\n--- Diagnostic breakdown ---\n")

for (tp in unique(all_pairs$Time)) {
  sub <- all_pairs[all_pairs$Time == tp, , drop = FALSE]

  cat("\nTimepoint:", tp, "\n")
  cat("Total genes tested:      ", nrow(sub), "\n")
  cat("Stroke_sig:              ", sum(sub$Stroke_sig, na.rm = TRUE), "\n")
  cat("NeuDep_sig:              ", sum(sub$NeuDep_sig, na.rm = TRUE), "\n")
  cat("Kept (union):            ", sum(sub$Keep, na.rm = TRUE), "\n")
}

cat("\nUnion signature size (unique genes kept across timepoints): ",
    length(unique(all_pairs$Geneid[all_pairs$Keep])), "\n")

cat("\n-----------------------------\n")

# -----------------------
# Summary
# -----------------------
cat("\nStroke inflammation signature (mouse) built.\n\n")

cat("Included timepoints: ",
    paste(c(if (use_3h) "3h" else NULL, if (use_24h) "24h" else NULL), collapse = ", "),
    "\n", sep = "")

cat("Stroke vs Sham threshold: ", use_stat, " <= ", p_stroke,
    " and |logFC| >= ", sprintf("%.3f", logfc_stroke), "\n", sep = "")

cat("NeuDep vs Stroke threshold: ", use_stat, " <= ", p_neudep,
    " and |logFC| >= ", sprintf("%.3f", logfc_neudep), "\n", sep = "")

cat("Inclusion rule: union of Stroke_sig and NeuDep_sig\n")

cat("\nMouse signature size: ", length(mouse_sig), "\n", sep = "")
cat("Mouse signature saved to:\n  ", mouse_out, "\n", sep = "")
cat("Signature evidence table saved to:\n  ", tsv_out, "\n", sep = "")

cat("\n")