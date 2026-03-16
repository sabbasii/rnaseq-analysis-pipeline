#!/usr/bin/env Rscript
# ============================================================
# scripts/dge/pde/04_pde_subgroup_consistency_summary.R
#
# Purpose
#   Within each subgroup (Age and Sex), run DIRECT contrasts:
#     Anti-Ly6G vs Stroke at:
#       - 3hPost-Reperfusion
#       - 24hPost-Stroke
#
#   Then extract PDE genes and produce a "consistency table":
#     Geneid | Young_3h | Aged_3h | Male_3h | Female_3h | Young_24h | ...
#   where each cell includes direction + significance (e.g., Down*, Up).
#
# Inputs
#   dge/inputs/matrices/counts_matrix.tsv
#   dge/inputs/metadata/metadata.tsv
#   metadata/gene_sets/pde_genes.txt
#
# Outputs
#   results/dge/pde/subgroups/
#     - pde_subgroup_consistency.long.tsv
#     - pde_subgroup_consistency.wide.tsv
#     - subgroup_edger_tables/<subgroup>_AntiLy6G_vs_Stroke_<tp>.all.tsv
#
# Significance definition
#   Use either --fdr or --pval (mutually exclusive) AND |logFC| >= log2(fc)
#   Defaults: --fdr 0.10 --fc 1.20
#
# Example
#   Rscript scripts/dge/pde/04_pde_subgroup_consistency_summary.R
#   `Rscript scripts/dge/pde/04_pde_subgroup_consistency_summary.R --fdr 0.05 --fc 1.2`
#   Rscript scripts/dge/pde/04_pde_subgroup_consistency_summary.R --pval 0.5  --fc 1.10
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
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

fdr_cut  <- get_arg("--fdr",  NA)
pval_cut <- get_arg("--pval", NA)
fc_cut   <- as.numeric(get_arg("--fc", "1.20"))

if (!is.na(fdr_cut))  fdr_cut  <- as.numeric(fdr_cut)
if (!is.na(pval_cut)) pval_cut <- as.numeric(pval_cut)

if (!is.na(fdr_cut) && !is.na(pval_cut)) stop("Use only one of --fdr or --pval (not both).")
if (is.na(fdr_cut) && is.na(pval_cut)) fdr_cut <- 0.10
if (is.na(fc_cut) || fc_cut <= 1) stop("--fc must be > 1")

logfc_thr <- log2(fc_cut)
use_stat <- if (!is.na(pval_cut)) "PValue" else "FDR"
stat_cut <- if (!is.na(pval_cut)) pval_cut else fdr_cut

# -----------------------
# Paths
# -----------------------
repo_root <- Sys.getenv("REPO_ROOT")
if (repo_root == "") repo_root <- normalizePath(getwd())

counts_file <- file.path(repo_root, "dge", "inputs", "matrices", "counts_matrix.tsv")
meta_file   <- file.path(repo_root, "dge", "inputs", "metadata", "metadata.tsv")
pde_file    <- file.path(repo_root, "metadata", "gene_sets", "pde_genes.txt")

need_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p)
need_file(counts_file)
need_file(meta_file)
need_file(pde_file)

out_dir <- file.path(repo_root, "results", "dge", "pde", "subgroups")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tab_out_dir <- file.path(out_dir, "subgroup_edger_tables")
dir.create(tab_out_dir, recursive = TRUE, showWarnings = FALSE)

out_long <- file.path(out_dir, "pde_subgroup_consistency.long.tsv")
out_wide <- file.path(out_dir, "pde_subgroup_consistency.wide.tsv")

# -----------------------
# Load inputs
# -----------------------
counts_dt <- fread(counts_file)
if (!("Geneid" %in% names(counts_dt))) {
  stop("counts_matrix.tsv must contain a 'Geneid' column.")
}
genes <- counts_dt$Geneid
sample_cols <- setdiff(names(counts_dt), "Geneid")

counts_mat <- as.matrix(counts_dt[, sample_cols, with = FALSE])
rownames(counts_mat) <- genes

meta <- fread(meta_file, data.table = FALSE)
req_cols <- c("SampleID", "treatment", "timePoint", "age", "sex")
miss <- setdiff(req_cols, colnames(meta))
if (length(miss) > 0) stop("metadata.tsv missing columns: ", paste(miss, collapse = ", "))

for (cc in req_cols) meta[[cc]] <- as.character(meta[[cc]])

# Keep only metadata samples that exist in counts
meta <- meta[meta$SampleID %in% colnames(counts_mat), , drop = FALSE]
if (nrow(meta) == 0) stop("No metadata SampleID values match counts matrix columns.")

# Reorder counts columns to match metadata order
counts_mat <- counts_mat[, meta$SampleID, drop = FALSE]

# -----------------------
# Load PDE genes
# -----------------------
pde <- fread(pde_file, header = FALSE, data.table = FALSE)[[1]]
pde <- trimws(as.character(pde))
pde <- pde[nzchar(pde)]
pde <- pde[!grepl("^#", pde)]
pde <- unique(pde)
if (length(pde) == 0) stop("No PDE genes loaded from: ", pde_file)

# -----------------------
# Helper: build simple group names + run edgeR
# -----------------------
run_subgroup_edger <- function(sub_name, meta_sub, counts_sub) {

  # Keep only the 4 groups we need for DIRECT contrasts
  meta_sub$group <- NA_character_
  meta_sub$group[meta_sub$treatment == "Anti-Ly6G" & meta_sub$timePoint == "3hPost-Reperfusion"] <- "anti_3h"
  meta_sub$group[meta_sub$treatment == "Stroke"    & meta_sub$timePoint == "3hPost-Reperfusion"] <- "stroke_3h"
  meta_sub$group[meta_sub$treatment == "Anti-Ly6G" & meta_sub$timePoint == "24hPost-Stroke"]     <- "anti_24h"
  meta_sub$group[meta_sub$treatment == "Stroke"    & meta_sub$timePoint == "24hPost-Stroke"]     <- "stroke_24h"

  keep_rows <- !is.na(meta_sub$group)
  meta_sub <- meta_sub[keep_rows, , drop = FALSE]
  if (nrow(meta_sub) == 0) {
    message("Skipping subgroup '", sub_name, "' (no Anti-Ly6G/Stroke samples at 3h/24h)")
    return(NULL)
  }

  # Subset counts to the same samples
  counts_sub <- counts_sub[, meta_sub$SampleID, drop = FALSE]

  meta_sub$group <- factor(meta_sub$group, levels = c("stroke_3h", "anti_3h", "stroke_24h", "anti_24h"))

  needed <- c("stroke_3h", "anti_3h", "stroke_24h", "anti_24h")
  if (!all(needed %in% levels(meta_sub$group))) {
    missing <- setdiff(needed, levels(meta_sub$group))
    message("Skipping subgroup '", sub_name, "' (missing groups: ", paste(missing, collapse = ", "), ")")
    return(NULL)
  }

  # edgeR
  y <- DGEList(counts = counts_sub)
  keep_genes <- filterByExpr(y, group = meta_sub$group)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  design <- model.matrix(~ 0 + group, data = meta_sub)
  colnames(design) <- gsub("^group", "", colnames(design))

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)

  c3  <- makeContrasts(anti_3h  - stroke_3h,  levels = design)
  c24 <- makeContrasts(anti_24h - stroke_24h, levels = design)

  tt3  <- topTags(glmQLFTest(fit, contrast = c3),  n = Inf, sort.by = "none")$table
  tt24 <- topTags(glmQLFTest(fit, contrast = c24), n = Inf, sort.by = "none")$table

  tt3$Geneid  <- rownames(tt3)
  tt24$Geneid <- rownames(tt24)

  keep_cols <- c("Geneid", "logFC", "PValue", "FDR")
  tt3  <- tt3[,  keep_cols, drop = FALSE]
  tt24 <- tt24[, keep_cols, drop = FALSE]

  f3  <- file.path(tab_out_dir, paste0(sub_name, "_AntiLy6G_vs_Stroke_3h.all.tsv"))
  f24 <- file.path(tab_out_dir, paste0(sub_name, "_AntiLy6G_vs_Stroke_24h.all.tsv"))
  fwrite(tt3,  f3,  sep = "\t")
  fwrite(tt24, f24, sep = "\t")

  list(tt3 = tt3, tt24 = tt24, files = c(f3, f24))
}

# -----------------------
# Subgroups
# -----------------------
subgroups <- list(
  Young  = function(m) m$age == "young",
  Aged   = function(m) m$age == "aged",
  Male   = function(m) m$sex == "male",
  Female = function(m) m$sex == "female"
)

# -----------------------
# Run
# -----------------------
all_rows <- list()

for (nm in names(subgroups)) {
  idx <- subgroups[[nm]](meta)
  meta_sub <- meta[idx, , drop = FALSE]
  if (nrow(meta_sub) < 6) {
    message("Skipping subgroup '", nm, "' (too few samples: ", nrow(meta_sub), ")")
    next
  }

  counts_sub <- counts_mat[, meta_sub$SampleID, drop = FALSE]

  res <- run_subgroup_edger(nm, meta_sub, counts_sub)
  if (is.null(res)) next

  for (tp in c("3h", "24h")) {
    tt <- if (tp == "3h") res$tt3 else res$tt24
    tt <- tt[tt$Geneid %in% pde, , drop = FALSE]
    if (nrow(tt) == 0) next

    sig_stat <- if (use_stat == "PValue") {
      !is.na(tt$PValue) & tt$PValue <= pval_cut
    } else {
      !is.na(tt$FDR) & tt$FDR <= fdr_cut
    }
    sig_fc <- !is.na(tt$logFC) & abs(tt$logFC) >= logfc_thr

    tt$Sig <- sig_stat & sig_fc
    tt$Dir <- ifelse(tt$logFC > 0, "Up", "Down")

    tt$Subgroup <- nm
    tt$TimePoint <- tp

    all_rows[[length(all_rows) + 1]] <- tt
  }
}

if (length(all_rows) == 0) {
  stop("No subgroup results were generated. Check subgroup sample availability at 3h/24h for Anti-Ly6G vs Stroke.")
}

long <- rbindlist(all_rows, fill = TRUE)
setcolorder(long, c("Geneid", "Subgroup", "TimePoint", "logFC", "PValue", "FDR", "Sig", "Dir"))
fwrite(long, out_long, sep = "\t")

# Wide consistency table
long$Cell <- paste0(long$Dir, ifelse(long$Sig, "*", ""))
long$key <- paste0(long$Subgroup, "_", long$TimePoint)

wide <- dcast(long, Geneid ~ key, value.var = "Cell", fill = "")

wide_lfc <- dcast(long, Geneid ~ key, value.var = "logFC", fill = NA_real_)
setnames(wide_lfc, old = names(wide_lfc)[-1], new = paste0(names(wide_lfc)[-1], "_logFC"))

wide_p <- dcast(long, Geneid ~ key, value.var = use_stat, fill = NA_real_)
setnames(wide_p, old = names(wide_p)[-1], new = paste0(names(wide_p)[-1], "_", use_stat))

wide_out <- merge(wide, wide_lfc, by = "Geneid", all = TRUE)
wide_out <- merge(wide_out, wide_p, by = "Geneid", all = TRUE)

fwrite(wide_out, out_wide, sep = "\t")

cat("Saved:\n  ", out_long, "\n  ", out_wide, "\n", sep = "")
cat("Cutoffs: ", use_stat, " <= ", stat_cut,
    " and |logFC| >= ", sprintf("%.3f", logfc_thr), "\n", sep = "")
cat("Subgroup edgeR tables written to:\n  ", tab_out_dir, "\n", sep = "")