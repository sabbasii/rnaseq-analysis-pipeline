#!/usr/bin/env Rscript
# ============================================================
# scripts/translational/build_mouse_to_human_ortholog_map_mgi.R
#
# Purpose
#   Build a mouse_symbol ↔ human_symbol ortholog mapping table from
#   MGI homology data. The resulting table is used by:
#
#     scripts/translational/map_mouse_signature_to_human.R
#
# Source
#   MGI homology report:
#     HOM_MouseHumanSequence.rpt
#
#   This tab-delimited report contains homology class identifiers and
#   gene symbols for mouse and human orthologs.
#
# Why MGI?
#   The Mouse Genome Informatics (MGI) database provides a widely used
#   and curated resource for mouse–human orthology relationships. Using
#   the downloadable report ensures the mapping step is reproducible.
#
# Inputs
#   Option A (default)
#     The script downloads the MGI homology report automatically.
#
#   Option B
#     Use a locally downloaded report:
#
#       --download false
#       --mgi_rpt <path>
#
# Outputs
#   metadata/orthologs/mouse_to_human_symbols.tsv
#
#   Columns:
#     class_key
#     mouse_symbol
#     human_symbol
#     mouse_mgi_id
#     human_hgnc_id
#
#   Each row represents one ortholog relationship. One-to-many or
#   many-to-one mappings are preserved as separate rows.
#
# Run
#
#   Default (download report automatically)
#
#     Rscript scripts/translational/build_mouse_to_human_ortholog_map_mgi.R
#
#   Using a local MGI report
#
#     Rscript scripts/translational/build_mouse_to_human_ortholog_map_mgi.R \
#       --download false \
#       --mgi_rpt metadata/orthologs/HOM_MouseHumanSequence.rpt
#
# Notes
#   - Orthology relationships are defined by shared homology class keys
#     in the MGI report.
#   - The script preserves all mouse–human pairs without collapsing
#     one-to-many relationships.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) stop(flag, " requires a value")
  args[i + 1]
}

repo_root <- Sys.getenv("REPO_ROOT")
if (!nzchar(repo_root)) repo_root <- normalizePath(getwd())

need_file <- function(path, label = path) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
}

orth_dir <- file.path(repo_root, "metadata", "orthologs")
dir.create(orth_dir, recursive = TRUE, showWarnings = FALSE)

out_map <- get_arg("--out", file.path(orth_dir, "mouse_to_human_symbols.tsv"))

mgi_url <- get_arg(
  "--mgi_url",
  "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
)

mgi_rpt <- get_arg("--mgi_rpt", file.path(orth_dir, "HOM_MouseHumanSequence.rpt"))

do_download <- tolower(get_arg("--download", "true")) %in% c("true","t","1","yes","y")

if (do_download) {
  message("Downloading MGI report:\n  ", mgi_url, "\n→ ", mgi_rpt)
  tryCatch({
    download.file(mgi_url, destfile = mgi_rpt, quiet = TRUE, mode = "wb")
  }, error = function(e) {
    stop("Download failed. Download manually and re-run with --download false --mgi_rpt <path>\n",
         "Error: ", e$message)
  })
}

need_file(mgi_rpt, "MGI homology report")

# Force tab parsing (this report is tab-delimited)
df <- fread(
  mgi_rpt,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  quote = "",
  fill = TRUE
)

cn0 <- colnames(df)
cn  <- tolower(cn0)

pick_col <- function(exacts = character(0), patterns = character(0)) {
  if (length(exacts) > 0) {
    ix <- which(cn %in% exacts)
    if (length(ix) > 0) return(ix[1])
  }
  if (length(patterns) > 0) {
    for (pat in patterns) {
      ix <- which(grepl(pat, cn, perl = TRUE))
      if (length(ix) > 0) return(ix[1])
    }
  }
  NA_integer_
}

ix_class <- pick_col(exacts = c("db class key", "db_class_key"), patterns = c("class.*key"))
ix_org   <- pick_col(exacts = c("common organism name", "common_organism_name"), patterns = c("organism"))
ix_sym   <- pick_col(exacts = c("symbol"), patterns = c("^symbol$"))
ix_mgi   <- pick_col(exacts = c("mouse mgi id", "mgi id", "mgi_id"), patterns = c("mgi.*id"))
ix_hgnc  <- pick_col(exacts = c("hgnc id", "hgnc_id"), patterns = c("hgnc.*id"))

if (is.na(ix_class) || is.na(ix_org) || is.na(ix_sym)) {
  stop("Could not find required columns in MGI report.\n",
       "Need: DB Class Key, Common Organism Name, Symbol\n",
       "Found columns: ", paste(cn0, collapse = ", "))
}

class_col <- cn0[ix_class]
org_col   <- cn0[ix_org]
sym_col   <- cn0[ix_sym]

mgi_col  <- if (!is.na(ix_mgi))  cn0[ix_mgi]  else NULL
hgnc_col <- if (!is.na(ix_hgnc)) cn0[ix_hgnc] else NULL

# Normalize
df[[class_col]] <- trimws(as.character(df[[class_col]]))
df[[org_col]]   <- trimws(as.character(df[[org_col]]))
df[[sym_col]]   <- trimws(as.character(df[[sym_col]]))

org_lc <- tolower(df[[org_col]])

# IMPORTANT: works for labels like "mouse, laboratory" and "human"
is_mouse <- grepl("^mouse", org_lc)
is_human <- grepl("^human", org_lc)

message("\n[MGI DEBUG] Organism examples (top 8):")
message(paste(utils::head(unique(df[[org_col]]), 8), collapse = " | "))
message("[MGI DEBUG] Mouse rows: ", sum(is_mouse, na.rm = TRUE),
        "  Human rows: ", sum(is_human, na.rm = TRUE))

df2 <- df[(is_mouse | is_human) & nzchar(df[[class_col]]) & nzchar(df[[sym_col]]), , drop = FALSE]

mouse_rows <- df2[grepl("^mouse", tolower(df2[[org_col]])), , drop = FALSE]
human_rows <- df2[grepl("^human", tolower(df2[[org_col]])), , drop = FALSE]

# Minimal tables
mouse_min <- mouse_rows[, c(class_col, sym_col, mgi_col), drop = FALSE]
human_min <- human_rows[, c(class_col, sym_col, hgnc_col), drop = FALSE]

colnames(mouse_min)[1:2] <- c("class_key", "mouse_symbol")
colnames(human_min)[1:2] <- c("class_key", "human_symbol")
if (!is.null(mgi_col))  colnames(mouse_min)[3] <- "mouse_mgi_id"
if (!is.null(hgnc_col)) colnames(human_min)[3] <- "human_hgnc_id"

# Pair within homology class (keeps one-to-many transparently)
pairs <- merge(mouse_min, human_min, by = "class_key", all = FALSE, sort = FALSE)

pairs$mouse_symbol <- trimws(as.character(pairs$mouse_symbol))
pairs$human_symbol <- trimws(as.character(pairs$human_symbol))
pairs <- pairs[nzchar(pairs$mouse_symbol) & nzchar(pairs$human_symbol), , drop = FALSE]
pairs <- unique(pairs)

fwrite(pairs, out_map, sep = "\t", quote = FALSE, na = "NA")

cat("\nMGI mouse→human ortholog mapping built.\n\n")
cat("MGI report used:\n  ", mgi_rpt, "\n", sep = "")
cat("Pairs written:\n  ", out_map, "\n\n", sep = "")
cat("Rows (mouse↔human pairs): ", nrow(pairs), "\n", sep = "")
cat("Unique mouse symbols:     ", length(unique(pairs$mouse_symbol)), "\n", sep = "")
cat("Unique human symbols:     ", length(unique(pairs$human_symbol)), "\n\n", sep = "")