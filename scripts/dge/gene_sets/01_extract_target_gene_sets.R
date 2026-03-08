#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

repo_root <- Sys.getenv("REPO_ROOT")
if (repo_root == "") repo_root <- getwd()

# Input
diff_file <- file.path(repo_root,
  "results/dge/model_treat_time/summary/diffGenes.edger.tsv")

# Gene lists
pde_file <- file.path(repo_root,
  "metadata/gene_sets/pde_pathway_genes.txt")

neut_file <- file.path(repo_root,
  "metadata/gene_sets/immune_neutrophil_genes.txt")

# Output
out_dir <- file.path(repo_root,
  "results/dge/target_gene_sets")

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Load data
diffGenes <- fread(diff_file)

pde <- fread(pde_file, header=FALSE)[[1]]
neut <- fread(neut_file, header=FALSE)[[1]]

# Extract
pde_res  <- diffGenes[Geneid %in% pde]
neut_res <- diffGenes[Geneid %in% neut]

# Save
fwrite(pde_res,
  file.path(out_dir, "PDE_genes_logFC.tsv"),
  sep="\t")

fwrite(neut_res,
  file.path(out_dir, "Neutrophil_genes_logFC.tsv"),
  sep="\t")

cat("Saved:\n", out_dir, "\n")

# Rscript scripts/dge/08_extract_target_gene_sets.R
