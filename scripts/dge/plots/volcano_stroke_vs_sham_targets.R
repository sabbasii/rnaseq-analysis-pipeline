#!/usr/bin/env Rscript
# ============================================================
# scripts/dge/plots/volcano_stroke_vs_sham_targets.R
#
# Purpose
#   Volcano plot for DIRECT Stroke vs Sham differential expression
#   at matched timepoints:
#     - stroke_vs_sham_3h
#     - stroke_vs_sham_24h
#
#   Highlights and labels target genes from:
#     - metadata/gene_sets/pde_genes.txt
#     - metadata/gene_sets/immune_neutrophil_genes.txt
#
# Inputs
#   - results/dge/model_direct_stroke_vs_sham/tables/<contrast>.all.tsv
#   - metadata/gene_sets/pde_genes.txt
#   - metadata/gene_sets/immune_neutrophil_genes.txt
#
# Output
#   - results/dge/target_gene_sets/volcano/volcano_<contrast>_targets.png
#
# Significance
#   Use either --fdr or --pval (mutually exclusive), plus fold-change:
#     statistic <= cutoff AND |logFC| >= log2(fc)
#
# Defaults (exploratory)
#   --contrast stroke_vs_sham_3h
#   --fdr 0.10
#   --fc  1.20
#   --label_mode sig_pde
#   --label_top 25
#   --sample_bg 0
#
# Label modes
#   top_targets : label top N significant targets (union of PDE + neutrophil targets)
#   sig_targets : label all significant targets (union)
#   sig_pde     : label all significant PDE genes only
#   all_pde     : label all PDE genes (significant or not)
#
# Examples
#   Rscript scripts/dge/plots/volcano_stroke_vs_sham_targets.R
#
#   Rscript scripts/dge/plots/volcano_stroke_vs_sham_targets.R \
#     --contrast stroke_vs_sham_24h --pval 0.05 --fc 1.25 --label_mode sig_targets
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

contrast   <- get_arg("--contrast", "stroke_vs_sham_3h")
fc_cut     <- as.numeric(get_arg("--fc", "1.20"))
label_mode <- get_arg("--label_mode", "sig_pde")
label_top  <- as.integer(get_arg("--label_top", "25"))
sample_bg  <- as.integer(get_arg("--sample_bg", "0"))

fdr_cut  <- get_arg("--fdr",  NA)
pval_cut <- get_arg("--pval", NA)
if (!is.na(fdr_cut))  fdr_cut  <- as.numeric(fdr_cut)
if (!is.na(pval_cut)) pval_cut <- as.numeric(pval_cut)

if (!is.na(fdr_cut) && !is.na(pval_cut)) stop("Use only one of --fdr or --pval (not both).")
if (is.na(fdr_cut) && is.na(pval_cut)) fdr_cut <- 0.10

if (is.na(fc_cut) || fc_cut <= 1) stop("--fc must be > 1 (e.g., 1.20)")
if (is.na(label_top) || label_top < 0) stop("--label_top must be >= 0")
if (is.na(sample_bg) || sample_bg < 0) stop("--sample_bg must be >= 0")

logfc_thr <- log2(fc_cut)

use_stat <- if (!is.na(pval_cut)) "PValue" else "FDR"
stat_cut <- if (!is.na(pval_cut)) pval_cut else fdr_cut

# -----------------------
# Paths
# -----------------------
repo_root <- Sys.getenv("REPO_ROOT")
if (repo_root == "") repo_root <- normalizePath(getwd())

tab_file <- file.path(
  repo_root,
  "results", "dge", "model_direct_stroke_vs_sham", "tables",
  paste0(contrast, ".all.tsv")
)

pde_file <- file.path(repo_root, "metadata", "gene_sets", "pde_genes.txt")
target_file <- file.path(repo_root, "metadata", "gene_sets", "immune_neutrophil_genes.txt")

out_dir <- file.path(repo_root, "results", "dge", "target_gene_sets", "volcano")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_png <- file.path(out_dir, paste0("volcano_", contrast, "_targets.png"))

need_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}
need_file(tab_file)
need_file(pde_file)
need_file(target_file)

# -----------------------
# Load gene sets
# -----------------------
read_genes <- function(path) {
  x <- fread(path, header = FALSE, data.table = FALSE)[[1]]
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  x <- x[!grepl("^#", x)]
  unique(x)
}

pde_genes <- read_genes(pde_file)
targets2  <- read_genes(target_file)
target_genes <- unique(c(pde_genes, targets2))

# -----------------------
# Load edgeR table
# -----------------------
df <- fread(tab_file, data.table = FALSE)
req <- c("Geneid", "logFC", "PValue", "FDR")
miss <- setdiff(req, colnames(df))
if (length(miss) > 0) stop("edgeR table missing columns: ", paste(miss, collapse = ", "))

df$Geneid <- as.character(df$Geneid)

# Y-axis
df$neglog10P   <- -log10(pmax(df$PValue, 1e-300))
df$neglog10FDR <- -log10(pmax(df$FDR, 1e-300))
df$Y <- if (use_stat == "PValue") df$neglog10P else df$neglog10FDR

# Significance
sig_stat <- if (use_stat == "PValue") {
  !is.na(df$PValue) & df$PValue <= pval_cut
} else {
  !is.na(df$FDR) & df$FDR <= fdr_cut
}
sig_fc <- !is.na(df$logFC) & abs(df$logFC) >= logfc_thr
sig <- sig_stat & sig_fc

# Regulation classes
df$Regulation <- "Not significant"
df$Regulation[sig & df$logFC >=  logfc_thr] <- "Upregulated"
df$Regulation[sig & df$logFC <= -logfc_thr] <- "Downregulated"
df$Regulation <- factor(df$Regulation, levels = c("Upregulated", "Downregulated", "Not significant"))

# Targets
df$is_target <- df$Geneid %in% target_genes
is_pde <- df$Geneid %in% pde_genes

# Optional background downsampling (always keep targets)
if (sample_bg > 0) {
  set.seed(1)
  bg <- df[!df$is_target, , drop = FALSE]
  if (nrow(bg) > sample_bg) bg <- bg[sample(seq_len(nrow(bg)), sample_bg), , drop = FALSE]
  df <- rbind(df[df$is_target, , drop = FALSE], bg)
}

# -----------------------
# Labels
# -----------------------
df$Label <- NA_character_

if (label_mode == "top_targets") {

  if (label_top > 0) {
    dft <- df[df$is_target & sig, , drop = FALSE]
    if (nrow(dft) > 0) {
      statv <- if (use_stat == "PValue") dft$PValue else dft$FDR
      dft$rankv <- statv
      dft <- dft[order(dft$rankv, -abs(dft$logFC)), , drop = FALSE]
      top_genes <- head(dft$Geneid, label_top)
      df$Label[df$Geneid %in% top_genes] <- df$Geneid[df$Geneid %in% top_genes]
    }
  }

} else if (label_mode == "sig_targets") {

  df$Label[df$is_target & sig] <- df$Geneid[df$is_target & sig]

} else if (label_mode == "sig_pde") {

  df$Label[is_pde & sig] <- df$Geneid[is_pde & sig]

} else if (label_mode == "all_pde") {

  df$Label[is_pde] <- df$Geneid[is_pde]

} else {
  stop("--label_mode must be one of: top_targets, sig_targets, sig_pde, all_pde")
}

# -----------------------
# Plot
# -----------------------
base_family <- "Helvetica"

p <- ggplot(df, aes(x = logFC, y = Y)) +

  # Main points
  geom_point(
    aes(color = Regulation),
    size = 0.7,
    alpha = 0.85
  ) +

  # Highlight significant target genes only
  geom_point(
    data = df[df$is_target & sig, , drop = FALSE],
    shape = 21,
    fill = "#55595b",
    color = "#55595b",
    stroke = 0.6,
    size = 1.0,
    inherit.aes = TRUE,
    show.legend = FALSE
  ) +

  # Threshold lines
  geom_vline(
    xintercept = c(-logfc_thr, logfc_thr),
    linetype = "dashed",
    linewidth = 0.4
  ) +

  geom_hline(
    yintercept = -log10(stat_cut),
    linetype = "dashed",
    linewidth = 0.4
  ) +

  # Labels with connector lines
  ggrepel::geom_text_repel(
    data = df[!is.na(df$Label), , drop = FALSE],
    aes(label = Label),
    size = 4.2,
    family = base_family,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.35,
    min.segment.length = 0,
    segment.color = "#55595b",
    segment.size = 0.4,
    segment.alpha = 0.9,
    seed = 1
  ) +

  # Custom colors (YOUR defined colors)
  scale_color_manual(
    values = c(
      "Upregulated"   = "#ed9696",
      "Downregulated" = "#a2bed3",
      "Not significant" = "grey75"
    )
  ) +

  # Labels and titles
  labs(
    title = contrast,
    subtitle = paste0(
      "Direct Stroke vs Sham | ",
      use_stat, " ≤ ", stat_cut,
      " and |logFC| ≥ log2(", fc_cut, ") = ",
      sprintf("%.3f", logfc_thr)
    ),
    x = "log2 fold change (logFC)",
    y = if (use_stat == "PValue") "-log10(PValue)" else "-log10(FDR)",
    color = "Regulation"
  ) +

  # Theme
  theme_bw(base_size = 12, base_family = base_family) +

  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey92")
  )

ggsave(out_png, p, width = 8, height = 6, dpi = 300)

cat("Saved:\n  ", out_png, "\n", sep = "")
cat("Input:\n  ", tab_file, "\n", sep = "")
cat("Targets (union of):\n  ", pde_file, "\n  ", target_file, "\n", sep = "")
cat("Cutoffs: ", use_stat, " <= ", stat_cut, " and |logFC| >= ", sprintf("%.3f", logfc_thr), "\n", sep = "")
cat("Label mode: ", label_mode, "\n", sep = "")


# Option 1 — Label top N targets (PDE + neutrophil) --> labels most significant targets only
# You can change number of labels

# Rscript scripts/dge/plots/volcano_stroke_vs_sham_targets.R \
#   --contrast stroke_vs_sham_3h \
#   --pval 0.05 \
#   --fc 1.2 \
#   --label_mode top_targets \
#   --label_top 25
# -----------------------------------------------
# Option 2 — Label ALL significant targets (PDE + neutrophil) --> Labels every significant PDE and neutrophil gene (recommended) **
# Rscript scripts/dge/plots/volcano_stroke_vs_sham_targets.R \
#   --contrast stroke_vs_sham_3h \
#   --pval 0.05 \
#   --fc 1.2 \
#   --label_mode sig_targets
# -----------------------------------------------
# Option 3 — Label only significant PDE genes --> Labels only PDE family genes
# Rscript scripts/dge/plots/volcano_stroke_vs_sham_targets.R \
#   --contrast stroke_vs_sham_3h \
#   --pval 0.05 \
#   --fc 1.2 \
#   --label_mode sig_pde
# -----------------------------------------------
# Option 4 — Label ALL PDE genes (significant or not)
# Rscript scripts/dge/plots/volcano_stroke_vs_sham_targets.R \
#   --contrast stroke_vs_sham_3h \
#   --pval 0.05 \
#   --fc 1.2 \
#   --label_mode all_pde
# -----------------------------------------------
# Optional modifier (can combine with any option) --> Limit background points for a cleaner plot)
# --sample_bg 5000