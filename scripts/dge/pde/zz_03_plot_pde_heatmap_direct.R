#!/usr/bin/env Rscript
# ============================================================
# scripts/dge/pde/03_plot_pde_heatmap_direct.R
#
# Purpose
#   Create a PDE-only heatmap for DIRECT NeuDep vs Stroke:
#     - neudep_vs_stroke_3h
#     - neudep_vs_stroke_24h
#
#   IMPORTANT UPDATE:
#     Only show PDE genes from pde_genes.txt that are FOUND in the
#     edgeR tables (i.e., present in expression/annotation used by edgeR).
#
#   Rows: PDE genes present in tables
#   Columns: 3h, 24h
#   Values: logFC
#   Significance marker: "*" on cells that pass thresholds
#
# Inputs
#   - results/dge/model_direct_neudep_vs_stroke/tables/neudep_vs_stroke_3h.all.tsv
#   - results/dge/model_direct_neudep_vs_stroke/tables/neudep_vs_stroke_24h.all.tsv
#   - metadata/gene_sets/pde_genes.txt
#
# Output
#   - results/dge/pde/heatmap/heatmap_pde_direct_neudep_vs_stroke.png
#
# Significance
#   Use either --fdr or --pval (mutually exclusive) AND
#   |logFC| >= log2(fc)
#
# Defaults (exploratory)
#   --fdr 0.10
#   --fc  1.20
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

tab_dir <- file.path(repo_root, "results", "dge", "model_direct_neudep_vs_stroke", "tables")
tab_3h  <- file.path(tab_dir, "neudep_vs_stroke_3h.all.tsv")
tab_24h <- file.path(tab_dir, "neudep_vs_stroke_24h.all.tsv")

pde_file <- file.path(repo_root, "metadata", "gene_sets", "pde_genes.txt")

out_dir <- file.path(repo_root, "results", "dge", "pde", "heatmap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_png <- file.path(out_dir, "heatmap_pde_direct_neudep_vs_stroke.png")

need_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}
need_file(tab_3h)
need_file(tab_24h)
need_file(pde_file)

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
# Read edgeR tables
# -----------------------
read_edger <- function(path, time_label) {
  df <- fread(path, data.table = FALSE)
  req <- c("Geneid", "logFC", "PValue", "FDR")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", path, ": ", paste(miss, collapse = ", "))

  df <- df[, req]
  df$Geneid <- as.character(df$Geneid)
  df$Time <- time_label
  df
}

df3  <- read_edger(tab_3h,  "3h")
df24 <- read_edger(tab_24h, "24h")
df_all <- rbind(df3, df24)

# -----------------------
# Keep ONLY PDE genes found in the tables (i.e., in expression data)
# -----------------------
genes_in_tables <- unique(df_all$Geneid)
pde_found <- intersect(pde, genes_in_tables)

if (length(pde_found) == 0) {
  msg <- paste0(
    "None of the PDE genes from pde_genes.txt were found in the edgeR tables.\n",
    "Check if pde_genes.txt uses a different ID type than Geneid (e.g., symbols vs Ensembl).\n",
    "Example Geneid values in tables:\n  ",
    paste(head(genes_in_tables, 10), collapse = ", "),
    "\n"
  )
  stop(msg)
}

df <- df_all[df_all$Geneid %in% pde_found, , drop = FALSE]

# Ensure we have exactly one row per (Geneid, Time)
dup_key <- duplicated(df[, c("Geneid", "Time")])
if (any(dup_key)) {
  # If edgeR table somehow has duplicates, keep the most significant row per key.
  # (Prefer smaller p-value / FDR depending on mode)
  stat_col <- if (use_stat == "PValue") "PValue" else "FDR"
  df <- as.data.table(df)
  setorder(df, Geneid, Time, get(stat_col), -abs(logFC))
  df <- unique(df, by = c("Geneid", "Time"))
  df <- as.data.frame(df)
}

# -----------------------
# Significance flag
# -----------------------
sig_stat <- if (use_stat == "PValue") {
  !is.na(df$PValue) & df$PValue <= pval_cut
} else {
  !is.na(df$FDR) & df$FDR <= fdr_cut
}
sig_fc <- !is.na(df$logFC) & abs(df$logFC) >= logfc_thr
df$Sig <- sig_stat & sig_fc
df$Star <- ifelse(df$Sig, "*", "")

# -----------------------
# Ordering (nice visual grouping)
# -----------------------
ord <- tapply(df$logFC, df$Geneid, function(x) max(abs(x), na.rm = TRUE))
ord[is.infinite(ord)] <- 0
gene_order <- names(sort(ord, decreasing = TRUE))

df$Geneid <- factor(df$Geneid, levels = gene_order)
df$Time <- factor(df$Time, levels = c("3h", "24h"))

# -----------------------
# Plot
# -----------------------
base_family <- "Helvetica"

p <- ggplot(df, aes(x = Time, y = Geneid, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = Star), size = 5, family = base_family) +
  labs(
    title = "PDE genes: Direct NeuDep vs Stroke",
    subtitle = paste0(
      "Cells marked * pass: ", use_stat, " ≤ ", stat_cut,
      " and |logFC| ≥ log2(", fc_cut, ") = ", sprintf("%.3f", logfc_thr)
    ),
    x = NULL,
    y = NULL,
    fill = "logFC"
  ) +
  theme_bw(base_size = 12, base_family = base_family) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 7, hjust = 0.5),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

ggsave(
  out_png, p,
  width = 3.5,
  height = max(4, 0.25 * length(pde_found) + 1.5),
  dpi = 300
)

cat("Saved:\n  ", out_png, "\n", sep = "")
cat("Included PDE genes found in tables: ", length(pde_found), " / ", length(pde), "\n", sep = "")
cat("Cutoffs: ", use_stat, " <= ", stat_cut, " and |logFC| >= ", sprintf("%.3f", logfc_thr), "\n", sep = "")

# Run
# Rscript scripts/dge/pde/03_plot_pde_heatmap_direct.R --pval 0.05 --fc 1.5