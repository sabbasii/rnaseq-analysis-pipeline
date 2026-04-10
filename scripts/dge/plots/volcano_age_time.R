#!/usr/bin/env Rscript
# ------------------------------------------------------------
# scripts/dge/plots/volcano_age_time.R
#
# Purpose
#   Volcano plots for age-by-time contrasts from:
#     results/dge/model_age_time/tables/*.all.tsv
#
# Thresholds (as requested)
#   - fcThreshold = log2(1.5)
#   - pvThreshold = 0.05
#   - label: top 25 up + top 25 down (by smallest PValue) beyond fcThreshold
#
# Inputs
#   - results/dge/model_age_time/tables/age_at_reperfusion.all.tsv
#   - results/dge/model_age_time/tables/age_at_followup.all.tsv
#
# Outputs
#   - results/dge/model_age_time/plots/volcano_age_3h.png
#   - results/dge/model_age_time/plots/volcano_age_24h.png
#
# Run (from repo root)
#   Rscript scripts/dge/plots/volcano_age_time.R
# ------------------------------------------------------------

# -------------------------------
# Install required packages if missing
# -------------------------------
required_pkgs <- c("ggrepel", "data.table", "ggplot2")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
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
    return(normalizePath(file.path(script_dir, "..", "..", "..")))
  }

  wd
}

need_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
}

make_volcano <- function(df, out_png) {
  req <- c("Geneid", "logFC", "PValue")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Missing columns in table: ", paste(miss, collapse = ", "))

  plotData <- df[df$PValue < 0.75, , drop = FALSE]

  fcThreshold <- log2(1.5)
  pvThreshold <- 0.05

  labelled_up <- plotData[plotData$logFC > fcThreshold, , drop = FALSE]
  if (nrow(labelled_up) > 0) labelled_up <- labelled_up[order(labelled_up$PValue), , drop = FALSE]
  labelled_up <- head(labelled_up, 25)

  labelled_dn <- plotData[plotData$logFC < -fcThreshold, , drop = FALSE]
  if (nrow(labelled_dn) > 0) labelled_dn <- labelled_dn[order(labelled_dn$PValue), , drop = FALSE]
  labelled_dn <- head(labelled_dn, 25)

  labelled <- rbind(labelled_up, labelled_dn)

  p <- ggplot(plotData, aes(x = logFC, y = -log10(PValue))) +
    geom_point(alpha = 0.5) +
    geom_text_repel(
      data = labelled,
      aes(label = Geneid),
      max.overlaps = Inf
    ) +
    geom_vline(xintercept = c(-fcThreshold, fcThreshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(pvThreshold), linetype = "dashed") +
    theme(text = element_text(size = 16))

  ggsave(filename = out_png, plot = p, width = 8, height = 6, dpi = 300)
}

# ---------- Paths ----------
repo_root <- resolve_repo_root()

tbl_dir <- file.path(repo_root, "results", "dge", "model_age_time", "tables")
plots_dir <- file.path(repo_root, "results", "dge", "model_age_time", "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

rep_path <- file.path(tbl_dir, "age_at_reperfusion.all.tsv")
fol_path <- file.path(tbl_dir, "age_at_followup.all.tsv")

need_file(rep_path, "age_at_reperfusion.all.tsv")
need_file(fol_path, "age_at_followup.all.tsv")

rep_df <- fread(rep_path, data.table = FALSE)
fol_df <- fread(fol_path, data.table = FALSE)

make_volcano(rep_df, file.path(plots_dir, "volcano_age_3h.png"))
make_volcano(fol_df, file.path(plots_dir, "volcano_age_24h.png"))

cat("Done. Wrote volcano plots to:\n  ", plots_dir, "\n", sep = "")