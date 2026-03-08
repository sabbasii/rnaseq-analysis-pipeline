#!/usr/bin/env Rscript
# ============================================================
# scripts/dge/pde/01_pde_neudep_vs_stroke_summary.R
#
# Purpose
#   Summarize PDE family differential expression for:
#     - NeuDep vs Stroke at 3h  (treatment_at_reperfusion)
#     - NeuDep vs Stroke at 24h (treatment_at_followup)
#
#   Pulls all PDE genes, keeps logFC/FDR, and flags "significant"
#   using a less strict cutoff (default):
#     FDR <= 0.10 and |logFC| >= log2(1.20)
#
# Inputs
#   - results/dge/model_treat_time/tables/treatment_at_reperfusion.all.tsv
#   - results/dge/model_treat_time/tables/treatment_at_followup.all.tsv
#   - metadata/gene_sets/pde_pathway_genes.txt  (or pde_genes.txt)
#
# Outputs
#   - results/dge/pde/neudep_vs_stroke_pde_summary.tsv
#
# Options
#   --fdr <value>   default 0.10
#   --fc  <value>   default 1.20   (fold-change threshold; uses log2(fc))
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------
# Args
# -----------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) stop(flag, " requires a value")
  args[i + 1]
}

fdr_cut <- as.numeric(get_arg("--fdr", "0.10"))
fc_cut  <- as.numeric(get_arg("--fc",  "1.20"))

if (is.na(fdr_cut) || fdr_cut <= 0 || fdr_cut >= 1) stop("--fdr must be between 0 and 1")
if (is.na(fc_cut)  || fc_cut <= 1) stop("--fc must be > 1 (e.g., 1.20)")

logfc_thr <- log2(fc_cut)

# -----------------------
# Paths
# -----------------------
repo_root <- Sys.getenv("REPO_ROOT")
if (repo_root == "") repo_root <- normalizePath(getwd())

tab_3h <- file.path(repo_root, "results", "dge", "model_treat_time", "tables",
                    "treatment_at_reperfusion.all.tsv")
tab_24h <- file.path(repo_root, "results", "dge", "model_treat_time", "tables",
                     "treatment_at_followup.all.tsv")

# Prefer expanded PDE pathway list if present, otherwise fall back to pde_genes.txt
pde_pathway_file <- file.path(repo_root, "metadata", "gene_sets", "pde_genes.txt")
pde_simple_file  <- file.path(repo_root, "metadata", "gene_sets", "pde_genes.txt")
pde_file <- if (file.exists(pde_pathway_file)) pde_pathway_file else pde_simple_file

out_dir <- file.path(repo_root, "results", "dge", "pde")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_tsv <- file.path(out_dir, "neudep_vs_stroke_pde_summary.tsv")

need_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

need_file(tab_3h)
need_file(tab_24h)
need_file(pde_file)

# -----------------------
# Load PDE list
# -----------------------
pde <- fread(pde_file, header = FALSE, data.table = FALSE)[[1]]
pde <- trimws(as.character(pde))
pde <- pde[nzchar(pde)]
pde <- pde[!grepl("^#", pde)]   # allow comment lines
pde <- unique(pde)

# -----------------------
# Load edgeR tables
# -----------------------
read_edger <- function(path) {
  df <- fread(path, data.table = FALSE)
  need <- c("Geneid", "logFC", "FDR")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop("edgeR table missing columns: ", paste(miss, collapse = ", "), "\nFile: ", path)
  df$Geneid <- as.character(df$Geneid)
  df
}

df3  <- read_edger(tab_3h)
df24 <- read_edger(tab_24h)

# -----------------------
# Keep PDE genes + flag significance
# -----------------------
keep_cols <- function(df, suffix) {
  sub <- df[df$Geneid %in% pde, c("Geneid", "logFC", "FDR"), drop = FALSE]
  colnames(sub) <- c("Geneid", paste0("logFC_", suffix), paste0("FDR_", suffix))
  sub
}

pde3  <- keep_cols(df3,  "3h")
pde24 <- keep_cols(df24, "24h")

out <- merge(pde3, pde24, by = "Geneid", all = TRUE, sort = FALSE)

flag_sig <- function(logFC, FDR) {
  !is.na(FDR) & !is.na(logFC) & (FDR <= fdr_cut) & (abs(logFC) >= logfc_thr)
}

out$Sig_3h  <- flag_sig(out$logFC_3h,  out$FDR_3h)
out$Sig_24h <- flag_sig(out$logFC_24h, out$FDR_24h)

# Add direction (helpful for quick scanning)
dir_from_logfc <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x > 0, "Up", ifelse(x < 0, "Down", "0")))
}
out$Dir_3h  <- dir_from_logfc(out$logFC_3h)
out$Dir_24h <- dir_from_logfc(out$logFC_24h)

# Order: show significant first
out$AnySig <- (out$Sig_3h %in% TRUE) | (out$Sig_24h %in% TRUE)
out <- out[order(out$AnySig, out$Sig_3h, out$Sig_24h, -abs(out$logFC_3h), -abs(out$logFC_24h), out$Geneid, na.last = TRUE), ]
out$AnySig <- NULL

# Write
fwrite(out, out_tsv, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote PDE NeuDep vs Stroke summary:\n  ", out_tsv, "\n", sep = "")
cat("Cutoffs used: FDR <= ", fdr_cut, " and |logFC| >= log2(", fc_cut, ") = ", round(logfc_thr, 3), "\n", sep = "")
cat("PDE list used:\n  ", pde_file, "\n", sep = "")

# ============================================================
# Run it:

# Default (less strict cutoffs):
# Rscript scripts/dge/pde/01_pde_neudep_vs_stroke_summary.R

# Less/more strict cutoffs:
# Rscript scripts/dge/pde/01_pde_neudep_vs_stroke_summary.R --fdr 0.15 --fc 1.15
# ============================================================
