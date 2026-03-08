#!/usr/bin/env Rscript
# ============================================================
# scripts/dge/pde/04_run_edger_neudep_vs_stroke_direct_and_interactions.R
#
# Purpose
#   Direct NeuDep (treated) vs Stroke differential expression at matched timepoints,
#   focused on PDE genes, and formal subgroup-dependence (interaction) tests.
#
#   This mirrors 03_run_edger_stroke_vs_sham_direct_and_interactions.R, but compares:
#       treated_3  vs stroke_3
#       treated_24 vs stroke_24
#
#   It also tests whether the NeuDep vs Stroke effect depends on:
#     - Age (young vs aged): Treatment × Age interaction
#     - Sex (male vs female): Treatment × Sex interaction
#
# Inputs
#   - dge/inputs/matrices/counts_matrix.tsv
#   - dge/inputs/metadata/metadata.tsv
#   - metadata/gene_sets/pde_genes.txt
#
# Outputs (consistent layout)
#   A) Full edgeR tables (all genes) for re-use elsewhere:
#      results/dge/model_direct_neudep_vs_stroke_timepoint/tables/
#         neudep_vs_stroke_3h.all.tsv
#         neudep_vs_stroke_24h.all.tsv
#         neudep_vs_stroke_3h.age_interaction.all.tsv
#         neudep_vs_stroke_24h.age_interaction.all.tsv
#         neudep_vs_stroke_3h.sex_interaction.all.tsv
#         neudep_vs_stroke_24h.sex_interaction.all.tsv
#         neudep_vs_stroke_3h.age_interaction_term.all.tsv
#         neudep_vs_stroke_24h.age_interaction_term.all.tsv
#         neudep_vs_stroke_3h.sex_interaction_term.all.tsv
#         neudep_vs_stroke_24h.sex_interaction_term.all.tsv
#
#   B) PDE-focused summaries:
#      results/dge/pde/neudep_vs_stroke_direct_timepoint/
#         pde_neudep_vs_stroke_direct_summary.tsv
#         pde_neudep_vs_stroke_age_interactions.tsv
#         pde_neudep_vs_stroke_sex_interactions.tsv
#         tables/
#            pde_neudep_vs_stroke_3h.tsv
#            pde_neudep_vs_stroke_24h.tsv
#
# Significance flagging (for PDE summary only)
#   Use either --fdr or --pval (mutually exclusive), plus --fc:
#     statistic <= cutoff AND |logFC| >= log2(fc)
#
# Defaults (exploratory)
#   --fdr 0.10
#   --fc  1.20
#
# Run
#   Rscript scripts/dge/pde/04_run_edger_neudep_vs_stroke_direct_and_interactions.R
#
#   Rscript scripts/dge/pde/04_run_edger_neudep_vs_stroke_direct_and_interactions.R \
#     --pval 0.05 --fc 1.20
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
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
# Repo root + paths
# -----------------------
resolve_repo_root <- function() {
  env_root <- Sys.getenv("REPO_ROOT")
  if (nzchar(env_root)) return(normalizePath(env_root))
  normalizePath(getwd())
}

need_file <- function(path, label = path) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
}

repo_root <- resolve_repo_root()

counts_in <- file.path(repo_root, "dge", "inputs", "matrices", "counts_matrix.tsv")
meta_in   <- file.path(repo_root, "dge", "inputs", "metadata", "metadata.tsv")
pde_file  <- file.path(repo_root, "metadata", "gene_sets", "pde_genes.txt")

need_file(counts_in, "counts_matrix.tsv")
need_file(meta_in,   "metadata.tsv")
need_file(pde_file,  "pde_genes.txt")

# Full-model outputs (all genes)
model_dir  <- file.path(repo_root, "results", "dge", "model_direct_neudep_vs_stroke_timepoint")
tables_dir <- file.path(model_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# PDE-focused outputs
pde_out_dir <- file.path(repo_root, "results", "dge", "pde", "neudep_vs_stroke_direct_timepoint")
pde_tab_dir <- file.path(pde_out_dir, "tables")
dir.create(pde_tab_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------
# Normalizers (must match your metadata values)
# -----------------------
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
  ifelse(grepl("^young$", x0), "Young",
  ifelse(grepl("^aged$|^old$", x0), "Aged", NA_character_))
}

normalize_sex <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  ifelse(grepl("^male$|^m$", x0), "Male",
  ifelse(grepl("^female$|^f$", x0), "Female", NA_character_))
}

# -----------------------
# Load inputs
# -----------------------
counts_df <- fread(counts_in, data.table = FALSE)
meta_df   <- fread(meta_in,   data.table = FALSE)

if (!"Geneid" %in% colnames(counts_df)) stop("counts_matrix.tsv must contain Geneid column.")
if (!"SampleID" %in% colnames(meta_df)) stop("metadata.tsv must contain SampleID column.")
if (!"treatment" %in% colnames(meta_df)) stop("metadata.tsv must contain treatment column.")
if (!"timePoint" %in% colnames(meta_df)) stop("metadata.tsv must contain timePoint column.")

# age/sex are required for interaction tests
if (!("age" %in% colnames(meta_df)) && !("Age" %in% colnames(meta_df))) {
  stop("metadata.tsv is missing an age column (expected 'age' or 'Age').")
}
if (!("sex" %in% colnames(meta_df)) && !("Sex" %in% colnames(meta_df))) {
  stop("metadata.tsv is missing a sex column (expected 'sex' or 'Sex').")
}

age_col <- if ("age" %in% colnames(meta_df)) "age" else "Age"
sex_col <- if ("sex" %in% colnames(meta_df)) "sex" else "Sex"

# Align metadata rows to counts sample order
sample_ids <- setdiff(colnames(counts_df), "Geneid")
meta_df$SampleID <- as.character(meta_df$SampleID)

if (!all(sample_ids %in% meta_df$SampleID)) {
  missing <- setdiff(sample_ids, meta_df$SampleID)
  stop("Counts contain SampleIDs not found in metadata:\n  ", paste(missing, collapse = "\n  "))
}
meta_df <- meta_df[match(sample_ids, meta_df$SampleID), , drop = FALSE]
stopifnot(identical(meta_df$SampleID, sample_ids))

# Normalize key fields
meta_df$Treatment <- normalize_treatment(meta_df$treatment)
meta_df$Time_hr   <- normalize_timepoint(meta_df$timePoint)
meta_df$AgeGroup  <- normalize_age(meta_df[[age_col]])
meta_df$SexGroup  <- normalize_sex(meta_df[[sex_col]])

if (anyNA(meta_df$Treatment)) stop("Unmapped treatment values present in metadata.")
if (anyNA(meta_df$Time_hr))   stop("Unmapped timePoint values present in metadata.")
if (anyNA(meta_df$AgeGroup))  stop("Unmapped age values present in metadata (expected Young/Aged).")
if (anyNA(meta_df$SexGroup))  stop("Unmapped sex values present in metadata (expected Male/Female).")

meta_df$AgeGroup <- factor(meta_df$AgeGroup, levels = c("Young", "Aged"))
meta_df$SexGroup <- factor(meta_df$SexGroup, levels = c("Male", "Female"))

# Count matrix
count_mat <- as.matrix(counts_df[, sample_ids, drop = FALSE])
rownames(count_mat) <- counts_df$Geneid
storage.mode(count_mat) <- "integer"

# PDE gene list
pde <- fread(pde_file, header = FALSE, data.table = FALSE)[[1]]
pde <- trimws(as.character(pde))
pde <- pde[nzchar(pde)]
pde <- pde[!grepl("^#", pde)]
pde <- unique(pde)

# -----------------------
# Helpers: fit edgeR for one timepoint subset
# -----------------------
fit_and_test <- function(time_hr, mode = c("main", "age_int", "sex_int")) {

  mode <- match.arg(mode)

  # subset to stroke+treated at timepoint
  keep_samples <- which(meta_df$Time_hr == time_hr & meta_df$Treatment %in% c("stroke", "treated"))

  if (length(keep_samples) < 4) {
    stop("Not enough samples at time ", time_hr, " for stroke+treated.")
  }

  sub_counts <- count_mat[, keep_samples, drop = FALSE]
  sub_meta   <- meta_df[keep_samples, , drop = FALSE]

  # Create factors with stable baseline (stroke as reference)
  sub_meta$Treatment2 <- factor(sub_meta$Treatment, levels = c("stroke", "treated"))

  y <- DGEList(counts = sub_counts)
  y <- calcNormFactors(y, method = "TMM")

  # Filtering within subset
  k <- filterByExpr(y, group = sub_meta$Treatment2)
  y <- y[k, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  if (mode == "main") {
    # adjusted for age+sex (no interaction)
    design <- model.matrix(~ Treatment2 + AgeGroup + SexGroup, data = sub_meta)
    coef_int <- NULL

  } else if (mode == "age_int") {
    design <- model.matrix(~ Treatment2 * AgeGroup + SexGroup, data = sub_meta)
    coef_int <- "Treatment2treated:AgeGroupAged"

  } else {
    design <- model.matrix(~ Treatment2 * SexGroup + AgeGroup, data = sub_meta)
    coef_int <- "Treatment2treated:SexGroupFemale"
  }

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  # main treated effect is always Treatment2treated in these parameterizations
  qlf_main <- glmQLFTest(fit, coef = "Treatment2treated")
  tab_main <- topTags(qlf_main, n = Inf, adjust.method = "BH")$table
  tab_main <- as.data.frame(tab_main)
  tab_main$Geneid <- rownames(tab_main)
  tab_main <- tab_main[, c("Geneid", setdiff(colnames(tab_main), "Geneid")), drop = FALSE]

  tab_int <- NULL
  if (!is.null(coef_int)) {
    if (!coef_int %in% colnames(design)) stop("Expected interaction coefficient not in design: ", coef_int)
    qlf_int <- glmQLFTest(fit, coef = coef_int)
    tab_int <- topTags(qlf_int, n = Inf, adjust.method = "BH")$table
    tab_int <- as.data.frame(tab_int)
    tab_int$Geneid <- rownames(tab_int)
    tab_int <- tab_int[, c("Geneid", setdiff(colnames(tab_int), "Geneid")), drop = FALSE]
  }

  list(main = tab_main, interaction = tab_int)
}

# -----------------------
# Run for 3h and 24h
# -----------------------
times <- c(3, 24)

# store PDE summaries
pde_main_list <- list()
pde_age_int_list <- list()
pde_sex_int_list <- list()

flag_sig <- function(df) {
  stat_ok <- if (use_stat == "PValue") {
    (!is.na(df$PValue) & df$PValue <= pval_cut)
  } else {
    (!is.na(df$FDR) & df$FDR <= fdr_cut)
  }
  fc_ok <- !is.na(df$logFC) & abs(df$logFC) >= logfc_thr
  stat_ok & fc_ok
}

for (t in times) {

  # ---- main effect (NeuDep vs Stroke, adjusted) ----
  res_main <- fit_and_test(t, "main")$main

  out_all_main <- file.path(tables_dir, paste0("neudep_vs_stroke_", t, "h.all.tsv"))
  fwrite(res_main, out_all_main, sep = "\t", quote = FALSE, na = "NA")

  # For direct NeuDep vs Stroke:
  #   logFC > 0 => higher in treated (NeuDep) than Stroke
  #   logFC < 0 => lower in treated (NeuDep) than Stroke
  pde_main <- res_main[res_main$Geneid %in% pde, c("Geneid","logFC","PValue","FDR"), drop = FALSE]
  pde_main$Time <- paste0(t, "h")
  pde_main$Sig <- flag_sig(pde_main)
  pde_main$Dir <- ifelse(is.na(pde_main$logFC), NA_character_,
                         ifelse(pde_main$logFC > 0, "Up_in_Treated",
                                ifelse(pde_main$logFC < 0, "Down_in_Treated", "0")))
  pde_main_list[[paste0(t, "h")]] <- pde_main

  fwrite(pde_main, file.path(pde_tab_dir, paste0("pde_neudep_vs_stroke_", t, "h.tsv")),
         sep = "\t", quote = FALSE, na = "NA")

  # ---- age interaction ----
  res_age <- fit_and_test(t, "age_int")
  out_all_age <- file.path(tables_dir, paste0("neudep_vs_stroke_", t, "h.age_interaction.all.tsv"))
  fwrite(res_age$main, out_all_age, sep = "\t", quote = FALSE, na = "NA")

  out_all_age_int <- file.path(tables_dir, paste0("neudep_vs_stroke_", t, "h.age_interaction_term.all.tsv"))
  fwrite(res_age$interaction, out_all_age_int, sep = "\t", quote = FALSE, na = "NA")

  pde_age_int <- res_age$interaction[res_age$interaction$Geneid %in% pde, c("Geneid","logFC","PValue","FDR"), drop = FALSE]
  pde_age_int$Time <- paste0(t, "h")
  pde_age_int$Interaction <- "Treatment:Age (Aged vs Young difference-in-effect)"
  pde_age_int$Sig <- flag_sig(pde_age_int)
  pde_age_int_list[[paste0(t, "h")]] <- pde_age_int

  # ---- sex interaction ----
  res_sex <- fit_and_test(t, "sex_int")
  out_all_sex <- file.path(tables_dir, paste0("neudep_vs_stroke_", t, "h.sex_interaction.all.tsv"))
  fwrite(res_sex$main, out_all_sex, sep = "\t", quote = FALSE, na = "NA")

  out_all_sex_int <- file.path(tables_dir, paste0("neudep_vs_stroke_", t, "h.sex_interaction_term.all.tsv"))
  fwrite(res_sex$interaction, out_all_sex_int, sep = "\t", quote = FALSE, na = "NA")

  pde_sex_int <- res_sex$interaction[res_sex$interaction$Geneid %in% pde, c("Geneid","logFC","PValue","FDR"), drop = FALSE]
  pde_sex_int$Time <- paste0(t, "h")
  pde_sex_int$Interaction <- "Treatment:Sex (Female vs Male difference-in-effect)"
  pde_sex_int$Sig <- flag_sig(pde_sex_int)
  pde_sex_int_list[[paste0(t, "h")]] <- pde_sex_int
}

# -----------------------
# Write PDE combined summaries
# -----------------------
pde_main_all <- rbindlist(pde_main_list, use.names = TRUE, fill = TRUE)
pde_main_all <- pde_main_all[order(pde_main_all$Geneid, pde_main_all$Time), ]

pde_main_out <- file.path(pde_out_dir, "pde_neudep_vs_stroke_direct_summary.tsv")
fwrite(pde_main_all, pde_main_out, sep = "\t", quote = FALSE, na = "NA")

pde_age_all <- rbindlist(pde_age_int_list, use.names = TRUE, fill = TRUE)
pde_age_all <- pde_age_all[order(pde_age_all$Geneid, pde_age_all$Time), ]
pde_age_out <- file.path(pde_out_dir, "pde_neudep_vs_stroke_age_interactions.tsv")
fwrite(pde_age_all, pde_age_out, sep = "\t", quote = FALSE, na = "NA")

pde_sex_all <- rbindlist(pde_sex_int_list, use.names = TRUE, fill = TRUE)
pde_sex_all <- pde_sex_all[order(pde_sex_all$Geneid, pde_sex_all$Time), ]
pde_sex_out <- file.path(pde_out_dir, "pde_neudep_vs_stroke_sex_interactions.tsv")
fwrite(pde_sex_all, pde_sex_out, sep = "\t", quote = FALSE, na = "NA")

cat(
  "Done.\n",
  "Full tables (all genes):\n  ", tables_dir, "\n",
  "PDE outputs:\n  ", pde_out_dir, "\n\n",
  "PDE significance used: ", use_stat, " <= ", stat_cut,
  " and |logFC| >= log2(", fc_cut, ") = ", sprintf("%.3f", logfc_thr), "\n",
  sep = ""
)

# To Run:
# Rscript scripts/dge/pde/04_run_edger_neudep_vs_stroke_direct_and_interactions.R --pval 0.05 --fc 1.2