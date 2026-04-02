#!/usr/bin/env Rscript
# ============================================================
# Script: dumbbell_stroke_vs_sham_and_neudep_vs_stroke.R
#
# Purpose
#   Generate a dumbbell plot showing paired differential expression
#   effects for phosphodiesterase (PDE) genes across matched timepoints.
#
#   Each gene is represented by two endpoints:
#
#     • Left circle  : Stroke vs Sham logFC
#     • Right circle : NeuDep vs Stroke logFC
#     • Connecting line shows magnitude and direction of shift
#
#   Separate panels are shown for each timepoint.
#
# Inputs
#   results/dge/pde/stroke_vs_sham_direct/
#     pde_stroke_vs_sham_direct_summary.tsv
#
#   results/dge/pde/neudep_vs_stroke_direct_timepoint/
#     pde_neudep_vs_stroke_direct_summary.tsv
#
# Output
#   results/dge/pde/plots/
#     pde_dumbbell_stroke_vs_sham__neudep_vs_stroke.png
#
# Requirements
#   R packages:
#     data.table
#     ggplot2 (>= 3.4.0)
#     forcats
#     scales
#
# Run
#   Rscript scripts/dge/plots/pde/dumbbell_stroke_vs_sham_and_neudep_vs_stroke.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# -----------------------
# Resolve paths
# -----------------------
repo_root <- Sys.getenv("REPO_ROOT")
if (!nzchar(repo_root)) repo_root <- normalizePath(getwd())

stroke_file <- file.path(
  repo_root,
  "results/dge/pde/stroke_vs_sham_direct/pde_stroke_vs_sham_direct_summary.tsv"
)

neudep_file <- file.path(
  repo_root,
  "results/dge/pde/neudep_vs_stroke_direct_timepoint/pde_neudep_vs_stroke_direct_summary.tsv"
)

out_dir <- file.path(repo_root, "results/dge/pde/plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(
  out_dir,
  "pde_dumbbell_stroke_vs_sham__neudep_vs_stroke.png"
)

need_file <- function(p) {
  if (!file.exists(p)) stop("Missing input file: ", p)
}

need_file(stroke_file)
need_file(neudep_file)

# -----------------------
# Load data
# -----------------------
stroke <- fread(stroke_file)
neudep <- fread(neudep_file)

required <- c("Geneid", "logFC", "Time", "Sig", "FDR", "PValue")

if (!all(required %in% names(stroke))) {
  stop("Stroke file missing columns: ",
       paste(setdiff(required, names(stroke)), collapse = ", "))
}

if (!all(required %in% names(neudep))) {
  stop("NeuDep file missing columns: ",
       paste(setdiff(required, names(neudep)), collapse = ", "))
}

setnames(
  stroke,
  c("logFC", "Sig", "FDR", "PValue"),
  c("logFC_stroke", "Sig_stroke", "FDR_stroke", "P_stroke")
)

setnames(
  neudep,
  c("logFC", "Sig", "FDR", "PValue"),
  c("logFC_neudep", "Sig_neudep", "FDR_neudep", "P_neudep")
)

df <- merge(
  stroke[, .(Geneid, Time, logFC_stroke)],
  neudep[, .(Geneid, Time, logFC_neudep)],
  by = c("Geneid", "Time"),
  all = FALSE
)

if (nrow(df) == 0) stop("No overlapping Geneid + Time entries.")

# Ensure consistent panel order
df$Time <- factor(df$Time, levels = c("3h", "24h"))

# Order genes by NeuDep effect (right endpoint)
df[, GeneOrder := forcats::fct_reorder(Geneid, logFC_neudep), by = Time]

# -----------------------
# Styling
# -----------------------
col_primary   <- "#a3c4dc"  # light blue (connector + left dot)
col_secondary <- "#0e668b"  # darker blue (right dot)

# -----------------------
# Plot
# -----------------------
p <- ggplot(df, aes(y = GeneOrder)) +

  geom_segment(
    aes(
      x = logFC_stroke,
      xend = logFC_neudep,
      yend = GeneOrder
    ),
    colour = col_primary,
    linewidth = 2.2
  ) +

  geom_segment(
    aes(
      x = -Inf,
      xend = logFC_stroke,
      yend = GeneOrder
    ),
    colour = "grey75",
    linewidth = 0.3
  ) +

  geom_point(
    aes(x = logFC_stroke),
    colour = col_primary,
    size = 3.6
  ) +

  geom_point(
    aes(x = logFC_neudep),
    colour = col_secondary,
    size = 3.6
  ) +
  geom_vline(
    xintercept = 0,
    colour = "red",
    linewidth = 0.6,
    linetype = "dashed"
  ) +

  facet_wrap(~Time, scales = "free_y") +

  labs(
    title = "PDE gene expression shifts across matched contrasts",
    subtitle = "Left: Stroke vs Sham   |   Right: NeuDep vs Stroke",
    x = "logFC",
    y = NULL
  ) +

  theme_minimal(base_size = 11) +

  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 11),

    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85"),

    axis.ticks = element_blank(),

    plot.margin = margin(6, 14, 6, 8),

    legend.position = "none"
  ) +

  scale_x_continuous(
    expand = expansion(mult = c(0.04, 0.06))
  )

# Reduced height for compact gene spacing
ggsave(
  out_png,
  p,
  width = 12,
  height = 3,
  dpi = 300
)

cat("Saved dumbbell plot:\n  ", out_png, "\n", sep = "")