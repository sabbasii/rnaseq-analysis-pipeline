#!/usr/bin/env Rscript
# ============================================================
# scripts/dge/plots/volcano_direct_neudep_vs_stroke_targets.R
#
# Purpose
#   Volcano plot for DIRECT NeuDep vs Stroke contrasts:
#     - neudep_vs_stroke_3h
#     - neudep_vs_stroke_24h
#
#   Highlights target genes from:
#     - metadata/gene_sets/pde_genes.txt
#     - metadata/gene_sets/target_pde_neutrophil.txt
#
# Inputs
#   - results/dge/model_direct_neudep_vs_stroke/tables/<contrast>.all.tsv
#   - metadata/gene_sets/pde_genes.txt
#   - metadata/gene_sets/target_pde_neutrophil.txt
#
# Output
#   - results/dge/target_gene_sets/volcano/volcano_<contrast>_targets.png
#
# Significance
#   Use either --fdr or --pval (mutually exclusive), plus fold-change:
#     statistic <= cutoff  AND  |logFC| >= log2(fc)
#
# Defaults
#   --contrast neudep_vs_stroke_3h
#   --fdr 0.10
#   --fc  1.20
#   --label_top 20
#   --label_mode sig_pde      # top_targets | sig_targets | sig_pde | all_pde
#   --sample_bg 0
#   --max_overlaps Inf        # ggrepel: set to Inf to allow many labels
#   --width  8
#   --height 6
#   --dpi 300
#
# Examples
#   Rscript scripts/dge/plots/volcano_direct_neudep_vs_stroke_targets.R
#
#   Rscript scripts/dge/plots/volcano_direct_neudep_vs_stroke_targets.R \
#     --contrast neudep_vs_stroke_24h --pval 0.05 --fc 1.25
#
#   # Label all significant PDEs at 3h:
#   Rscript scripts/dge/plots/volcano_direct_neudep_vs_stroke_targets.R \
#     --contrast neudep_vs_stroke_3h --label_mode sig_pde
#
#   # Label every PDE member no matter what:
#   Rscript scripts/dge/plots/volcano_direct_neudep_vs_stroke_targets.R \
#     --contrast neudep_vs_stroke_3h --label_mode all_pde
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
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

contrast   <- get_arg("--contrast", "neudep_vs_stroke_3h")
fc_cut     <- as.numeric(get_arg("--fc", "1.20"))
label_top  <- as.integer(get_arg("--label_top", "20"))
label_mode <- get_arg("--label_mode", "sig_pde")  # top_targets | sig_targets | sig_pde | all_pde
sample_bg  <- as.integer(get_arg("--sample_bg", "0"))

# ggrepel / output controls
max_overlaps <- get_arg("--max_overlaps", "Inf")
max_overlaps <- if (tolower(max_overlaps) == "inf") Inf else as.integer(max_overlaps)

out_w   <- as.numeric(get_arg("--width", "8"))
out_h   <- as.numeric(get_arg("--height", "6"))
out_dpi <- as.integer(get_arg("--dpi", "300"))

fdr_cut  <- get_arg("--fdr",  NA)
pval_cut <- get_arg("--pval", NA)
if (!is.na(fdr_cut))  fdr_cut  <- as.numeric(fdr_cut)
if (!is.na(pval_cut)) pval_cut <- as.numeric(pval_cut)

if (!is.na(fdr_cut) && !is.na(pval_cut)) stop("Use only one of --fdr or --pval (not both).")
if (is.na(fdr_cut) && is.na(pval_cut)) fdr_cut <- 0.10

if (is.na(fc_cut) || fc_cut <= 1) stop("--fc must be > 1 (e.g., 1.20)")
if (is.na(label_top) || label_top < 0) stop("--label_top must be >= 0")
if (is.na(sample_bg) || sample_bg < 0) stop("--sample_bg must be >= 0")
if (is.na(out_w) || out_w <= 0) stop("--width must be > 0")
if (is.na(out_h) || out_h <= 0) stop("--height must be > 0")
if (is.na(out_dpi) || out_dpi <= 0) stop("--dpi must be > 0")

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
  "results", "dge", "model_direct_neudep_vs_stroke", "tables",
  paste0(contrast, ".all.tsv")
)

pde_file <- file.path(repo_root, "metadata", "gene_sets", "pde_genes.txt")
target_file <- file.path(repo_root, "metadata", "gene_sets", "target_pde_neutrophil.txt")

out_dir <- file.path(repo_root, "results", "dge", "target_gene_sets", "volcano")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}
need_file(tab_file)
need_file(pde_file)
need_file(target_file)

# -----------------------
# Load targets
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
df$neglog10P <- -log10(pmax(df$PValue, 1e-300))
df$neglog10FDR <- -log10(pmax(df$FDR, 1e-300))
df$Y <- if (use_stat == "PValue") df$neglog10P else df$neglog10FDR

# Regulation category for color
df$Regulation <- "Not significant"

sig_stat <- if (use_stat == "PValue") {
  !is.na(df$PValue) & df$PValue <= pval_cut
} else {
  !is.na(df$FDR) & df$FDR <= fdr_cut
}

sig_fc <- !is.na(df$logFC) & abs(df$logFC) >= logfc_thr
sig <- sig_stat & sig_fc

df$Regulation[sig & df$logFC >=  logfc_thr] <- "Upregulated"
df$Regulation[sig & df$logFC <= -logfc_thr] <- "Downregulated"
df$Regulation <- factor(df$Regulation, levels = c("Upregulated", "Downregulated", "Not significant"))

df$is_target <- df$Geneid %in% target_genes

# Optional background downsampling (keep all targets)
if (sample_bg > 0) {
  set.seed(1)
  bg <- df[!df$is_target, , drop = FALSE]
  if (nrow(bg) > sample_bg) bg <- bg[sample(seq_len(nrow(bg)), sample_bg), , drop = FALSE]
  df <- rbind(df[df$is_target, , drop = FALSE], bg)
}

# -----------------------
# Labels (supports --label_mode)
# -----------------------
df$Label <- NA_character_

is_pde <- df$Geneid %in% pde_genes

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
# Plot (clean, standard)
# -----------------------
base_family <- "Helvetica"

p <- ggplot(df, aes(x = logFC, y = Y)) +
  geom_point(aes(color = Regulation), size = 0.7, alpha = 0.85) +  # dot size
  geom_point(
    data = df[df$is_target & sig, , drop = FALSE],  # Highlight only significant target genes
    shape = 21, fill = NA, color = "black", stroke = 0.6, size = 1.0, # outline dot size
    inherit.aes = TRUE, show.legend = FALSE
  ) +
  geom_vline(xintercept = c(-logfc_thr, logfc_thr), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(stat_cut), linetype = "dashed", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = df[!is.na(df$Label), , drop = FALSE],
    aes(label = Label),
    size = 3.5,
    family = base_family,
    fontface = "bold",
    max.overlaps = max_overlaps,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    seed = 1
  ) +
  # Change dot color
  scale_color_manual(
    values = c(
      "Upregulated" = "#ed9696",
      "Downregulated" = "#a2bed3",
      "Not significant" = "grey75"
    )
  ) +
  labs(
    title = contrast,
    subtitle = paste0("Direct NeuDep vs Stroke | ", use_stat, " ≤ ", stat_cut,
                      " and |logFC| ≥ log2(", fc_cut, ")"),
    x = "log2 fold change (logFC)",
    y = if (use_stat == "PValue") "-log10(PValue)" else "-log10(FDR)",
    color = "Regulation"
  ) +
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

out_png <- file.path(out_dir, paste0("volcano_", contrast, "_targets.png"))
ggsave(out_png, p, width = out_w, height = out_h, dpi = out_dpi)

cat("Saved:\n  ", out_png, "\n", sep = "")
cat("Targets highlighted (union of):\n  ", pde_file, "\n  ", target_file, "\n", sep = "")
cat("Cutoffs: ", use_stat, " <= ", stat_cut, " and |logFC| >= ", round(logfc_thr, 3), "\n", sep = "")
cat("Label mode: ", label_mode, "\n", sep = "")
cat("ggrepel max.overlaps: ", if (is.infinite(max_overlaps)) "Inf" else as.character(max_overlaps), "\n", sep = "")
cat("Output: ", out_w, "x", out_h, " inches @ ", out_dpi, " dpi\n", sep = "")


# Example Run:
# Rscript scripts/dge/plots/volcano_direct_neudep_vs_stroke_targets.R \
#     --contrast neudep_vs_stroke_3h --pval 0.05 --fc 1.25 --label_mode sig_targets