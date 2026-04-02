#!/usr/bin/env Rscript
# ------------------------------------------------------------
# scripts/dge/plots/volcano_treat_time.R
#
# Purpose
#   Volcano plots for treat_time contrasts (same style as old script).
#
# Inputs (from script 02)
#   results/dge/model_treat_time/tables/<contrast>.all.tsv
#
# Outputs
#   results/dge/model_treat_time/plots/volcano/<contrast>.png
#
# Run (from repo root)
#   Rscript scripts/dge/plots/volcano_treat_time.R
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
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
    return(normalizePath(file.path(script_dir, "..", "..", "..")))
  }
  wd
}

need_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
}

repo_root <- resolve_repo_root()

model_root <- file.path(repo_root, "results", "dge", "model_treat_time")
tables_dir <- file.path(model_root, "tables")
plots_dir  <- file.path(model_root, "plots", "volcano")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

contrasts <- c(
  "stroke_at_reperfusion",
  "stroke_at_followup",
  "treatment_at_reperfusion",
  "treatment_at_followup"
)

paths <- file.path(tables_dir, paste0(contrasts, ".all.tsv"))
for (i in seq_along(paths)) need_file(paths[i], paste0("table for ", contrasts[i]))

fcThreshold <- log2(1.5)
pvThreshold <- 0.05

plot_one <- function(contrast_name, in_path) {
  df <- fread(in_path, data.table = FALSE)

  needed <- c("Geneid", "logFC", "PValue", "FDR")
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) stop(basename(in_path), " is missing columns: ", paste(miss, collapse = ", "))

  plotData <- df[df$PValue < 0.75, , drop = FALSE]

  labelled <- plotData[
    (plotData$logFC < -2 & plotData$FDR <= 0.1) |
      (plotData$logFC > fcThreshold & plotData$FDR <= 0.1),
    , drop = FALSE
  ]

  p <- ggplot(plotData, aes(x = logFC, y = -log10(PValue))) +
    geom_point(alpha = 0.5) +
    geom_text_repel(
      data = labelled,
      aes(label = Geneid)
    ) +
    geom_vline(
      xintercept = c(-fcThreshold, fcThreshold),
      col = "gray",
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(pvThreshold),
      col = "red",
      linetype = "dashed"
    ) +
    labs() +
    scale_y_log10() +
    theme(text = element_text(size = 16))

  out_path <- file.path(plots_dir, paste0(contrast_name, ".png"))
  ggsave(filename = out_path, plot = p, width = 8, height = 6, dpi = 300)
  out_path
}

outs <- mapply(plot_one, contrasts, paths, USE.NAMES = FALSE)
cat("Wrote volcano plots to:\n", paste0("  ", outs, collapse = "\n"), "\n")
