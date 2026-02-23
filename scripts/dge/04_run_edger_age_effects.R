#!/usr/bin/env Rscript
# ------------------------------------------------------------
# scripts/dge/03_run_edger_age_effects.R
#
# Purpose
#   Run edgeR QL differential expression for AGE-stratified effects
#   using grouped factor: age_treat_time = age_Treatment_Time
#
# Inputs (from 01_prepare_edger_inputs.R)
#   - dge/inputs/matrices/counts_matrix.tsv
#   - dge/inputs/metadata/metadata.tsv
#
# Output
#   - results/dge/model_age_treat_time/tables/<contrast>.all.tsv
#   - results/dge/model_age_treat_time/tables_filtered/<contrast>.fdr25_logfc1p25.tsv
#   - results/dge/model_age_treat_time/summary/ageDiffGenes.tsv
#   - results/dge/model_age_treat_time/gene_lists/*.txt
#
# Run (from repo root)
#   Rscript scripts/dge/03_run_edger_age_effects.R
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
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
    return(normalizePath(file.path(script_dir, "..", "..")))
  }

  wd
}

need_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
}

normalize_treatment <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  ifelse(grepl("^sham$", x0), "sham",
  ifelse(grepl("^stroke$", x0) | grepl("^mcao$", x0), "stroke",
  ifelse(grepl("anti", x0) | grepl("ly6g", x0) | grepl("neude", x0) | grepl("^treated$", x0), "treated",
         NA_character_)))
}

normalize_timepoint <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  ifelse(grepl("baseline", x0) | grepl("^bl$", x0) | grepl("^0$", x0), 0,
  ifelse(grepl("3h", x0) | grepl("reperfusion", x0) | grepl("3hpr", x0), 3,
  ifelse(grepl("24h", x0) | grepl("follow", x0) | grepl("post.?stroke", x0), 24,
         NA_integer_)))
}

normalize_age <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  ifelse(grepl("^young$", x0) | grepl("^y$", x0), "young",
  ifelse(grepl("^aged$", x0) | grepl("^old$", x0) | grepl("^a$", x0), "aged",
         NA_character_))
}

repo_root <- resolve_repo_root()

counts_in <- file.path(repo_root, "dge", "inputs", "matrices", "counts_matrix.tsv")
meta_in   <- file.path(repo_root, "dge", "inputs", "metadata", "metadata.tsv")

need_file(counts_in, "counts_matrix.tsv")
need_file(meta_in,   "metadata.tsv")

model_dir   <- file.path(repo_root, "results", "dge", "model_age_treat_time")
tables_dir  <- file.path(model_dir, "tables")
filt_dir    <- file.path(model_dir, "tables_filtered")
sum_dir     <- file.path(model_dir, "summary")
lists_dir   <- file.path(model_dir, "gene_lists")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filt_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(lists_dir,  recursive = TRUE, showWarnings = FALSE)

counts_df <- fread(counts_in, data.table = FALSE)
meta_df   <- fread(meta_in,   data.table = FALSE)

if (!"Geneid" %in% colnames(counts_df)) stop("counts_matrix.tsv must contain Geneid column.")
if (!"SampleID" %in% colnames(meta_df)) stop("metadata.tsv must contain SampleID column.")
if (!"age" %in% colnames(meta_df)) stop("metadata.tsv must contain age column.")
if (!"treatment" %in% colnames(meta_df)) stop("metadata.tsv must contain treatment column.")
if (!"timePoint" %in% colnames(meta_df)) stop("metadata.tsv must contain timePoint column.")

sample_ids <- setdiff(colnames(counts_df), "Geneid")

meta_df$SampleID <- as.character(meta_df$SampleID)
if (!all(sample_ids %in% meta_df$SampleID)) {
  missing <- setdiff(sample_ids, meta_df$SampleID)
  stop("Counts contain SampleIDs not found in metadata:\n  ", paste(missing, collapse = "\n  "))
}
meta_df <- meta_df[match(sample_ids, meta_df$SampleID), , drop = FALSE]
stopifnot(identical(meta_df$SampleID, sample_ids))

meta_df$Treatment <- normalize_treatment(meta_df$treatment)
meta_df$Time_hr   <- normalize_timepoint(meta_df$timePoint)
meta_df$Age2      <- normalize_age(meta_df$age)

if (anyNA(meta_df$Treatment)) {
  bad <- unique(meta_df$treatment[is.na(meta_df$Treatment)])
  stop("Unmapped treatment values in metadata:\n  ", paste(bad, collapse = "\n  "))
}
if (anyNA(meta_df$Time_hr)) {
  bad <- unique(meta_df$timePoint[is.na(meta_df$Time_hr)])
  stop("Unmapped timePoint values in metadata:\n  ", paste(bad, collapse = "\n  "))
}
if (anyNA(meta_df$Age2)) {
  bad <- unique(meta_df$age[is.na(meta_df$Age2)])
  stop("Unmapped age values in metadata:\n  ", paste(bad, collapse = "\n  "))
}

meta_df$age_treat_time <- factor(paste(meta_df$Age2, meta_df$Treatment, meta_df$Time_hr, sep = "_"))

count_mat <- as.matrix(counts_df[, sample_ids, drop = FALSE])
rownames(count_mat) <- counts_df$Geneid
storage.mode(count_mat) <- "integer"

y <- DGEList(counts = count_mat, genes = data.frame(Geneid = counts_df$Geneid))

keep <- filterByExpr(y, group = meta_df$age_treat_time)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~0 + age_treat_time, data = meta_df)
colnames(design) <- levels(meta_df$age_treat_time)

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

needed_levels <- c(
  "young_sham_0","young_sham_3","young_sham_24",
  "young_stroke_0","young_stroke_3","young_stroke_24",
  "young_treated_0","young_treated_3","young_treated_24",
  "aged_sham_0","aged_sham_3","aged_sham_24",
  "aged_stroke_0","aged_stroke_3","aged_stroke_24",
  "aged_treated_0","aged_treated_3","aged_treated_24"
)

missing_lvls <- setdiff(needed_levels, colnames(design))
if (length(missing_lvls) > 0) {
  stop(
    "Design is missing required age_treat_time levels:\n  ",
    paste(missing_lvls, collapse = "\n  "),
    "\nCheck that metadata includes both ages across sham/stroke/treated at 0,3,24."
  )
}

contrast_mat <- makeContrasts(
  stroke_at_reperfusion_age_effect =
    (aged_stroke_3  - aged_stroke_0)  - (young_stroke_3 - young_stroke_0),

  stroke_at_followup_age_effect =
    (aged_stroke_24 - aged_stroke_3)  - (young_stroke_24 - young_stroke_3),

  treatment_at_reperfusion_age_effect =
    (aged_treated_3 - aged_treated_0) - (young_treated_3 - young_treated_0),

  stroke_at_followup_aged =
    (aged_stroke_24 - aged_stroke_3)  - (aged_sham_24 - aged_sham_3),

  treatment_at_followup_aged =
    (aged_treated_24 - aged_treated_3) - (aged_stroke_24 - aged_stroke_3),

  levels = design
)

LOGFC_CUTOFF <- log2(1.25)
FDR_CUTOFF <- 0.25

all_tables <- list()

run_one <- function(name) {
  qlf <- glmQLFTest(fit, contrast = contrast_mat[, name])
  tt <- topTags(qlf, n = Inf, adjust.method = "BH")$table
  tt <- as.data.frame(tt)
  tt$Geneid <- rownames(tt)
  tt <- tt[, c("Geneid", setdiff(colnames(tt), "Geneid")), drop = FALSE]

  out_all <- file.path(tables_dir, paste0(name, ".all.tsv"))
  fwrite(tt, out_all, sep = "\t", quote = FALSE, na = "NA")

  if (!all(c("FDR","logFC") %in% colnames(tt))) stop("Expected FDR and logFC columns in edgeR output.")
  tt_filt <- tt[!is.na(tt$FDR) & tt$FDR <= FDR_CUTOFF & !is.na(tt$logFC) & abs(tt$logFC) >= LOGFC_CUTOFF, , drop = FALSE]
  out_filt <- file.path(filt_dir, paste0(name, ".fdr25_logfc1p25.tsv"))
  fwrite(tt_filt, out_filt, sep = "\t", quote = FALSE, na = "NA")

  all_tables[[name]] <<- tt
}

invisible(lapply(colnames(contrast_mat), run_one))

ageDiffGenes <- data.frame(Geneid = sort(unique(unlist(lapply(all_tables, function(d) d$Geneid)))), stringsAsFactors = FALSE)
for (nm in names(all_tables)) {
  tmp <- all_tables[[nm]][, c("Geneid","logFC")]
  colnames(tmp) <- c("Geneid", nm)
  ageDiffGenes <- merge(ageDiffGenes, tmp, by = "Geneid", all.x = TRUE, sort = FALSE)
}
ageDiffGenes <- ageDiffGenes[rowSums(is.na(ageDiffGenes[, setdiff(colnames(ageDiffGenes), "Geneid"), drop = FALSE])) < (ncol(ageDiffGenes) - 1), , drop = FALSE]

fwrite(ageDiffGenes, file.path(sum_dir, "ageDiffGenes.tsv"), sep = "\t", quote = FALSE, na = "NA")

write_gene_list <- function(vec, filename) {
  vec <- unique(vec[!is.na(vec) & nzchar(vec)])
  fwrite(data.frame(Geneid = vec), file.path(lists_dir, filename), sep = "\t", quote = FALSE, col.names = FALSE)
}

if (all(c("stroke_at_reperfusion_age_effect","treatment_at_reperfusion_age_effect") %in% colnames(ageDiffGenes))) {
  restoredFunction <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_reperfusion_age_effect > 0 &
                                           ageDiffGenes$stroke_at_reperfusion_age_effect < 0]
  upregulated_treated_3h <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_reperfusion_age_effect > 0 &
                                                 is.na(ageDiffGenes$stroke_at_reperfusion_age_effect)]
  change_in_regulation_treated_3h <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_reperfusion_age_effect < 0 &
                                                          ageDiffGenes$stroke_at_reperfusion_age_effect > 0]
  downregulated_treated_3h <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_reperfusion_age_effect < 0 &
                                                   is.na(ageDiffGenes$stroke_at_reperfusion_age_effect)]
  furtherdown_treated_3h <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_reperfusion_age_effect <
                                                 ageDiffGenes$stroke_at_reperfusion_age_effect &
                                                 ageDiffGenes$stroke_at_reperfusion_age_effect < 0]

  write_gene_list(restoredFunction, "restoredFunction.txt")
  write_gene_list(upregulated_treated_3h, "upregulated_treated_3h.txt")
  write_gene_list(change_in_regulation_treated_3h, "change_in_regulation_treated_3h.txt")
  write_gene_list(downregulated_treated_3h, "downregulated_treated_3h.txt")
  write_gene_list(furtherdown_treated_3h, "furtherdown_treated_3h.txt")
}

if (all(c("stroke_at_followup_aged","treatment_at_followup_aged") %in% colnames(ageDiffGenes))) {
  restoredFunction_followup <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_followup_aged > 0 &
                                                    ageDiffGenes$stroke_at_followup_aged < 0]
  upregulated_treated_followup <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_followup_aged > 0 &
                                                        is.na(ageDiffGenes$stroke_at_followup_aged)]
  change_in_regulation_treated_followup <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_followup_aged < 0 &
                                                                 ageDiffGenes$stroke_at_followup_aged > 0]
  downregulated_treated_followup <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_followup_aged < 0 &
                                                          is.na(ageDiffGenes$stroke_at_followup_aged)]
  furtherdown_treated_followup <- ageDiffGenes$Geneid[ageDiffGenes$treatment_at_followup_aged <
                                                        ageDiffGenes$stroke_at_followup_aged &
                                                        ageDiffGenes$stroke_at_followup_aged < 0]

  write_gene_list(restoredFunction_followup, "restoredFunction_followup.txt")
  write_gene_list(upregulated_treated_followup, "upregulated_treated_followup.txt")
  write_gene_list(change_in_regulation_treated_followup, "change_in_regulation_treated_followup.txt")
  write_gene_list(downregulated_treated_followup, "downregulated_treated_followup.txt")
  write_gene_list(furtherdown_treated_followup, "furtherdown_treated_followup.txt")
}

cat("Done. Wrote outputs to:\n  ", model_dir, "\n", sep = "")
