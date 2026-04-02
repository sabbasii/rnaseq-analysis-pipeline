#!/usr/bin/env Rscript
# ============================================================
# scripts/translational/map_mouse_signature_to_human.R
#
# Purpose
#   Map a mouse gene signature list (one gene per line) to human symbols
#   using a user-provided ortholog mapping table (TSV/CSV).
#
#   This script filters the mapping table to ONLY the genes present in:
#     metadata/gene_sets/stroke_inflammation_signature_mouse.txt
#
# Inputs
#   1) Mouse signature (default):
#      metadata/gene_sets/stroke_inflammation_signature_mouse.txt
#
#   2) Ortholog mapping file (required):
#      --orthologs <path>
#      Must contain at least two columns:
#        mouse_symbol, human_symbol
#      (case-insensitive; common aliases accepted)
#
# Outputs
#   1) Human signature list:
#      metadata/gene_sets/stroke_inflammation_signature_human.txt
#
#   2) QC table (mapping evidence):
#      results/dge/target_gene_sets/signatures/stroke_inflammation_signature_mouse_to_human_map.tsv
#
# Run
  # Rscript scripts/translational/map_mouse_signature_to_human.R \
  #   --orthologs metadata/orthologs/mouse_to_human_symbols.tsv
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

# -----------------------
# Repo root + paths
# -----------------------
repo_root <- Sys.getenv("REPO_ROOT")
if (!nzchar(repo_root)) repo_root <- normalizePath(getwd())

need_file <- function(path, label = path) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
}

mouse_sig_path <- get_arg(
  "--mouse_sig",
  file.path(repo_root, "metadata", "gene_sets", "stroke_inflammation_signature_mouse.txt")
)

orthologs_path <- get_arg("--orthologs", NA)
if (is.na(orthologs_path)) stop("Required: --orthologs <mapping.tsv|csv>")

need_file(mouse_sig_path, "mouse signature .txt")
need_file(orthologs_path, "ortholog mapping file")

human_sig_path <- get_arg(
  "--human_out",
  file.path(repo_root, "metadata", "gene_sets", "stroke_inflammation_signature_human.txt")
)

qc_dir <- file.path(repo_root, "results", "dge", "target_gene_sets", "signatures")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

qc_out <- get_arg(
  "--qc_out",
  file.path(qc_dir, "stroke_inflammation_signature_mouse_to_human_map.tsv")
)

dir.create(dirname(human_sig_path), recursive = TRUE, showWarnings = FALSE)

# -----------------------
# Read mouse signature list
# -----------------------
mouse_sig <- readLines(mouse_sig_path, warn = FALSE)
mouse_sig <- trimws(mouse_sig)
mouse_sig <- mouse_sig[nzchar(mouse_sig)]
mouse_sig <- unique(mouse_sig)

if (length(mouse_sig) == 0) stop("Mouse signature list is empty: ", mouse_sig_path)

# -----------------------
# Read ortholog mapping file
# -----------------------
mapdf <- fread(orthologs_path, data.table = FALSE)
if (nrow(mapdf) == 0) stop("Ortholog mapping file is empty: ", orthologs_path)

cn0 <- colnames(mapdf)
cn <- tolower(cn0)

pick_col <- function(cands) {
  ix <- which(cn %in% cands)
  if (length(ix) == 0) return(NA_integer_)
  ix[1]
}

mouse_col <- pick_col(c("mouse_symbol","mouse","mgi_symbol","symbol_mouse","mousegene","gene_mouse"))
human_col <- pick_col(c("human_symbol","human","hgnc_symbol","symbol_human","humangene","gene_human"))

if (is.na(mouse_col) || is.na(human_col)) {
  stop(
    "Ortholog file must contain columns for mouse and human symbols.\n",
    "Accepted examples: mouse_symbol + human_symbol (recommended).\n",
    "Found columns: ", paste(cn0, collapse = ", ")
  )
}

mcol <- cn0[mouse_col]
hcol <- cn0[human_col]

mapdf[[mcol]] <- trimws(as.character(mapdf[[mcol]]))
mapdf[[hcol]] <- trimws(as.character(mapdf[[hcol]]))

mapdf <- mapdf[nzchar(mapdf[[mcol]]) & nzchar(mapdf[[hcol]]), , drop = FALSE]

# -----------------------
# FILTER mapping to only signature genes
# -----------------------
submap <- mapdf[mapdf[[mcol]] %in% mouse_sig, , drop = FALSE]

# Build QC table (keep all pairs for transparency)
qc <- data.frame(
  mouse_symbol = submap[[mcol]],
  human_symbol = submap[[hcol]],
  stringsAsFactors = FALSE
)

if (nrow(qc) > 0) {
  n_human_per_mouse <- aggregate(human_symbol ~ mouse_symbol, qc, function(x) length(unique(x)))
  colnames(n_human_per_mouse)[2] <- "n_human_symbols"
  qc <- merge(qc, n_human_per_mouse, by = "mouse_symbol", all.x = TRUE, sort = FALSE)
} else {
  qc$n_human_symbols <- integer(0)
}

mapped_mouse <- unique(qc$mouse_symbol)
unmapped_mouse <- setdiff(mouse_sig, mapped_mouse)

human_sig <- sort(unique(qc$human_symbol))
human_sig <- human_sig[nzchar(human_sig)]

# -----------------------
# Write outputs
# -----------------------
writeLines(human_sig, human_sig_path)

qc_unmapped <- data.frame(
  mouse_symbol = unmapped_mouse,
  human_symbol = NA_character_,
  n_human_symbols = 0L,
  stringsAsFactors = FALSE
)

qc_all <- rbind(qc[, c("mouse_symbol","human_symbol","n_human_symbols")], qc_unmapped)
qc_all <- qc_all[order(qc_all$mouse_symbol), , drop = FALSE]

fwrite(qc_all, qc_out, sep = "\t", quote = FALSE, na = "NA")

# -----------------------
# Summary
# -----------------------
cat("\nMouse→Human signature mapping complete.\n\n")
cat("Mouse signature genes:        ", length(mouse_sig), "\n", sep = "")
cat("Mouse genes mapped:           ", length(mapped_mouse), "\n", sep = "")
cat("Mouse genes unmapped:         ", length(unmapped_mouse), "\n", sep = "")
cat("Human signature (unique):     ", length(human_sig), "\n\n", sep = "")

cat("Mouse input:\n  ", mouse_sig_path, "\n", sep = "")
cat("Ortholog mapping used:\n  ", orthologs_path, "\n", sep = "")
cat("Human signature written:\n  ", human_sig_path, "\n", sep = "")
cat("QC mapping table written:\n  ", qc_out, "\n\n", sep = "")