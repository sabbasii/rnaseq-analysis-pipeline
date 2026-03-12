#!/usr/bin/env Rscript
# ============================================================
# scripts/dge/pde/05_pde_subgroup_consistency_report_and_interactions.R
#
# Outputs (TSV only)
#   results/dge/pde/subgroups/reports/
#     - pde_subgroup_consistency_summary.tsv
#     - pde_interaction_tests.tsv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

repo_root <- Sys.getenv("REPO_ROOT")
if (repo_root == "") repo_root <- normalizePath(getwd())

need_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p)

wide_file   <- file.path(repo_root, "results", "dge", "pde", "subgroups", "pde_subgroup_consistency.wide.tsv")
counts_file <- file.path(repo_root, "dge", "inputs", "matrices", "counts_matrix.tsv")
meta_file   <- file.path(repo_root, "dge", "inputs", "metadata", "metadata.tsv")
pde_file    <- file.path(repo_root, "metadata", "gene_sets", "pde_genes.txt")

need_file(wide_file); need_file(counts_file); need_file(meta_file); need_file(pde_file)

out_dir <- file.path(repo_root, "results", "dge", "pde", "subgroups", "reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_tsv <- file.path(out_dir, "pde_subgroup_consistency_summary.tsv")
out_int <- file.path(out_dir, "pde_interaction_tests.tsv")

# ------------------------------------------------------------
# 1) Consistency summary from wide table
# ------------------------------------------------------------
wide <- fread(wide_file)

dir_cols <- c("Aged_24h","Aged_3h","Female_24h","Female_3h","Male_24h","Male_3h","Young_24h","Young_3h")
dir_cols <- dir_cols[dir_cols %in% names(wide)]
if (!("Geneid" %in% names(wide))) stop("wide table missing Geneid")
if (length(dir_cols) == 0) stop("No direction columns found in wide table")

cell_clean <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "\"\"", "NA")] <- NA_character_
  x
}
get_dir <- function(x) {
  x <- cell_clean(x)
  ifelse(is.na(x), NA_character_, sub("\\*$", "", x))
}
get_sig <- function(x) {
  x <- cell_clean(x)
  !is.na(x) & grepl("\\*$", x)
}

summ <- wide[, {
  d <- unlist(lapply(.SD, get_dir), use.names = FALSE)
  s <- unlist(lapply(.SD, get_sig), use.names = FALSE)

  n_present <- sum(!is.na(d))
  n_up      <- sum(d == "Up", na.rm = TRUE)
  n_down    <- sum(d == "Down", na.rm = TRUE)
  n_sig     <- sum(s, na.rm = TRUE)

  direction_consistency <- if (n_present == 0) "NoData"
  else if (n_up == n_present) "AllUp"
  else if (n_down == n_present) "AllDown"
  else "Mixed"

  list(
    n_present = n_present,
    n_up = n_up,
    n_down = n_down,
    n_sig = n_sig,
    direction_consistency = direction_consistency,
    sig_fraction = ifelse(n_present > 0, n_sig / n_present, NA_real_)
  )
}, by = Geneid, .SDcols = dir_cols]

summ[, rank_score := fifelse(direction_consistency %in% c("AllDown","AllUp"),
                             1000 + 10*n_sig + n_present,
                             10*n_sig + n_present)]
setorder(summ, -rank_score, Geneid)
summ[, rank_score := NULL]

fwrite(summ, out_tsv, sep = "\t")

# ------------------------------------------------------------
# 2) Formal interaction tests (edgeR QL) per timepoint
#    Model: ~ treatment * covar
#    Test:  interaction coefficient(s) (overall F-test if multiple)
#
# NOTE: We *do not* enforce hard sample-size guards anymore.
#       Instead, we:
#         - try to fit
#         - auto-handle rank-deficient designs
#         - skip only when interaction is not estimable or fit/test fails
# ------------------------------------------------------------

# PDE list
pde <- fread(pde_file, header = FALSE, data.table = FALSE)[[1]]
pde <- trimws(as.character(pde))
pde <- pde[nzchar(pde)]
pde <- pde[!grepl("^#", pde)]
pde <- unique(pde)

# Counts
counts_dt <- fread(counts_file)
if (!("Geneid" %in% names(counts_dt))) stop("counts_matrix.tsv must have Geneid")
genes <- counts_dt$Geneid
sample_cols <- setdiff(names(counts_dt), "Geneid")
counts_mat <- as.matrix(counts_dt[, sample_cols, with = FALSE])
rownames(counts_mat) <- genes

# Metadata
meta <- fread(meta_file, data.table = FALSE)
req_cols <- c("SampleID","treatment","timePoint","age","sex")
miss <- setdiff(req_cols, colnames(meta))
if (length(miss) > 0) stop("metadata.tsv missing: ", paste(miss, collapse = ", "))

for (cc in req_cols) meta[[cc]] <- as.character(meta[[cc]])

# Match meta to counts
meta <- meta[meta$SampleID %in% colnames(counts_mat), , drop = FALSE]
if (nrow(meta) == 0) stop("No metadata SampleID match counts columns.")
counts_mat <- counts_mat[, meta$SampleID, drop = FALSE]

# Restrict to Stroke vs Anti-Ly6G (NeuDep) only
meta <- meta[meta$treatment %in% c("Stroke","Anti-Ly6G"), , drop = FALSE]
counts_mat <- counts_mat[, meta$SampleID, drop = FALSE]

# Utility: reduce a design matrix to full rank (drops redundant columns)
reduce_to_full_rank <- function(design) {
  q <- qr(design)
  r <- q$rank
  if (r < ncol(design)) {
    keep <- q$pivot[seq_len(r)]
    design <- design[, keep, drop = FALSE]
  }
  design
}

# Safe interaction runner
run_interaction <- function(tp, covar_name) {
  m <- meta[meta$timePoint == tp, , drop = FALSE]
  if (nrow(m) == 0) {
    message("SKIP ", tp, " ", covar_name, ": no samples")
    return(NULL)
  }

  # Subset counts
  cm <- counts_mat[, m$SampleID, drop = FALSE]

  # Factors (force treatment baseline)
  m$treatment <- factor(m$treatment, levels = c("Stroke","Anti-Ly6G"))
  m$covar <- factor(m[[covar_name]])

  # If only one level exists, interaction is not defined
  if (nlevels(droplevels(m$treatment)) < 2) {
    message("SKIP ", tp, " ", covar_name, ": 1 treatment level")
    return(NULL)
  }
  if (nlevels(droplevels(m$covar)) < 2) {
    message("SKIP ", tp, " ", covar_name, ": 1 covariate level")
    return(NULL)
  }

  design <- model.matrix(~ treatment * covar, data = m)
  design <- reduce_to_full_rank(design)
  cn <- colnames(design)

  # Identify treatment coef (main effect) if present
  trt_coef <- cn[grepl("^treatment", cn)]
  trt_coef <- trt_coef[grepl("Anti", trt_coef)]
  trt_coef <- if (length(trt_coef) >= 1) trt_coef[1] else NA_character_

  # Identify interaction term(s) if estimable
  int_terms <- cn[grepl(":", cn) & grepl("^treatment", cn)]
  int_terms <- int_terms[grepl("Anti", int_terms)]
  if (length(int_terms) < 1) {
    message("SKIP ", tp, " ", covar_name, ": no interaction term estimable (empty cell / rank drop)")
    return(NULL)
  }

  # Fit edgeR robustly
  fit <- tryCatch({
    y <- DGEList(counts = cm)
    keep <- filterByExpr(y, design)
    if (sum(keep) < 10) stop("too few expressed genes after filterByExpr")
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    glmQLFit(y, design, robust = TRUE)
  }, error = function(e) {
    message("SKIP ", tp, " ", covar_name, ": fit failed (", conditionMessage(e), ")")
    return(NULL)
  })
  if (is.null(fit)) return(NULL)

  # Interaction test (overall F-test if multiple coefs)
  int_idx <- match(int_terms, colnames(design))
  int_idx <- int_idx[!is.na(int_idx)]
  if (length(int_idx) < 1) {
    message("SKIP ", tp, " ", covar_name, ": interaction coef index missing after rank reduction")
    return(NULL)
  }

  tt_int <- tryCatch({
    topTags(glmQLFTest(fit, coef = int_idx), n = Inf, sort.by = "none")$table
  }, error = function(e) {
    message("SKIP ", tp, " ", covar_name, ": interaction test failed (", conditionMessage(e), ")")
    return(NULL)
  })
  if (is.null(tt_int)) return(NULL)
  tt_int$Geneid <- rownames(tt_int)

  # Treatment main effect test (if present)
  tt_trt <- NULL
  if (!is.na(trt_coef)) {
    trt_idx <- match(trt_coef, colnames(design))
    if (!is.na(trt_idx)) {
      tt_trt <- tryCatch({
        topTags(glmQLFTest(fit, coef = trt_idx), n = Inf, sort.by = "none")$table
      }, error = function(e) NULL)
      if (!is.null(tt_trt)) {
        tt_trt$Geneid <- rownames(tt_trt)
        tt_trt <- tt_trt[, c("Geneid","logFC","PValue","FDR"), drop = FALSE]
      }
    }
  }

  # Keep only essential interaction cols
  tt_int <- tt_int[, c("Geneid","logFC","PValue","FDR"), drop = FALSE]
  setnames(tt_int, c("logFC","PValue","FDR"), c("logFC_int","PValue_int","FDR_int"))

  # Merge treatment + interaction (treatment may be missing if dropped)
  if (is.null(tt_trt)) {
    out <- tt_int
    out$logFC <- NA_real_
    out$PValue <- NA_real_
    out$FDR <- NA_real_
    out <- out[, c("Geneid","logFC","PValue","FDR","logFC_int","PValue_int","FDR_int"), drop = FALSE]
  } else {
    out <- merge(tt_trt, tt_int, by = "Geneid", all = FALSE)
  }

  # Restrict to PDE
  out <- out[out$Geneid %in% pde, , drop = FALSE]
  if (nrow(out) == 0) {
    message("NOTE ", tp, " ", covar_name, ": model ran but no PDE genes present after filtering")
    return(NULL)
  }

  # Add metadata about what was tested
  out$timePoint <- tp
  out$covariate <- covar_name
  out$coef_treatment <- trt_coef
  out$coef_interaction <- paste(int_terms, collapse = ";")
  out$n_samples <- nrow(m)

  # Add cell table string so you can see sparsity quickly
  tab <- table(m$treatment, m$covar)
  out$cell_counts <- paste(capture.output(print(tab)), collapse = " | ")

  # Reorder columns
  out <- out[, c("Geneid","timePoint","covariate","n_samples",
                 "logFC","PValue","FDR",
                 "logFC_int","PValue_int","FDR_int",
                 "coef_treatment","coef_interaction",
                 "cell_counts"), drop = FALSE]
  out
}

# Choose timepoints of interest (your labels are correct)
timepoints <- c("3hPost-Reperfusion","24hPost-Stroke")

all_int <- list()
for (tp in timepoints) {
  a <- run_interaction(tp, "age")
  s <- run_interaction(tp, "sex")
  if (!is.null(a) && nrow(a) > 0) all_int[[length(all_int)+1]] <- as.data.table(a)
  if (!is.null(s) && nrow(s) > 0) all_int[[length(all_int)+1]] <- as.data.table(s)
}

# Always write a TSV WITH HEADERS (even if no rows)
template <- data.table(
  Geneid = character(),
  timePoint = character(),
  covariate = character(),
  n_samples = integer(),
  logFC = numeric(),
  PValue = numeric(),
  FDR = numeric(),
  logFC_int = numeric(),
  PValue_int = numeric(),
  FDR_int = numeric(),
  coef_treatment = character(),
  coef_interaction = character(),
  cell_counts = character()
)

if (length(all_int) == 0) {
  fwrite(template, out_int, sep = "\t")
} else {
  int_dt <- rbindlist(all_int, fill = TRUE)
  fwrite(int_dt, out_int, sep = "\t")
}

cat("Saved:\n", sep = "")
cat("  ", out_tsv, "\n", sep = "")
cat("  ", out_int, "\n", sep = "")

# Run it
# Rscript scripts/dge/pde/05_pde_subgroup_consistency_report_and_interactions.R