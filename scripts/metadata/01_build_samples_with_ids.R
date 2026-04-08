#!/usr/bin/env Rscript
# ------------------------------------------------------------
# scripts/metadata/01_build_samples_with_ids.R
#
# Purpose
#   Build standardized metadata table with:
#     - sampleID  (mouseID_treatment_timePoint)
#     - ageSex    (e.g., agedFemale, youngMale)
#     - age       (young/aged)
#     - sex       (male/female)
#
# Inputs
#   metadata/samples.tsv                        (NO HEADER; 4 columns:
#                                               position, mouseID, treatment, timePoint)
#   metadata/UpdatedSampleList_SimaProject.txt  (HAS HEADER; includes position + ageSex)
#
# Output
#   metadata/samples.with_ids.tsv
#
# Run (from repo root)
#   Rscript scripts/metadata/01_build_samples_with_ids.R
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

resolve_repo_root <- function() {
  env_root <- Sys.getenv("REPO_ROOT")
  if (nzchar(env_root)) return(normalizePath(env_root))

  wd <- normalizePath(getwd())
  if (file.exists(file.path(wd, "README.md")) || dir.exists(file.path(wd, "metadata"))) return(wd)

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    script_dir <- dirname(normalizePath(script_path))
    return(normalizePath(file.path(script_dir, "..", "..")))
  }
  wd
}

need_file <- function(path, label = path) {
  if (!file.exists(path)) stop("Missing file: ", label, "\n  ", path)
}

parse_age_sex <- function(ageSex) {
  x <- tolower(trimws(as.character(ageSex)))

  age <- ifelse(grepl("^young", x), "young",
         ifelse(grepl("^aged",  x), "aged", NA_character_))

  sex <- ifelse(grepl("female$", x), "female",
         ifelse(grepl("male$",   x), "male", NA_character_))

  data.frame(age = age, sex = sex, stringsAsFactors = FALSE)
}

build_sample_id <- function(mouseID, treatment, timePoint) {
  paste(trimws(as.character(mouseID)),
        trimws(as.character(treatment)),
        trimws(as.character(timePoint)),
        sep = "_")
}

repo_root <- resolve_repo_root()

samples_path <- file.path(repo_root, "metadata", "samples.tsv")
extra_path   <- file.path(repo_root, "metadata", "UpdatedSampleList_SimaProject.txt")
out_path     <- file.path(repo_root, "metadata", "samples.with_ids.tsv")

need_file(samples_path, "metadata/samples.tsv")
need_file(extra_path,   "metadata/UpdatedSampleList_SimaProject.txt")

samples <- fread(samples_path, header = FALSE, sep = "\t", data.table = FALSE)
if (ncol(samples) < 4) {
  stop("metadata/samples.tsv must have 4 tab-separated columns (no header).\nFound: ", ncol(samples))
}
samples <- samples[, 1:4, drop = FALSE]
colnames(samples) <- c("position", "mouseID", "treatment", "timePoint")

extra <- fread(extra_path, data.table = FALSE)

need_cols <- c("position", "ageSex")
missing_cols <- setdiff(need_cols, colnames(extra))
if (length(missing_cols) > 0) {
  stop(
    "UpdatedSampleList_SimaProject.txt is missing required columns:\n  ",
    paste(missing_cols, collapse = "\n  "),
    "\nColumns found:\n  ",
    paste(colnames(extra), collapse = "\n  ")
  )
}

extra_keep <- extra[, c("position", "ageSex"), drop = FALSE]

merged <- merge(
  samples,
  extra_keep,
  by = "position",
  all.x = TRUE,
  sort = FALSE
)

parsed <- parse_age_sex(merged$ageSex)
merged$age <- parsed$age
merged$sex <- parsed$sex

merged$sampleID <- build_sample_id(merged$mouseID, merged$treatment, merged$timePoint)

cat("\n========== QC REPORT (metadata build) ==========\n")

missing_join <- merged[is.na(merged$ageSex), , drop = FALSE]
if (nrow(missing_join) > 0) {
  cat("\nRows with NO ageSex match (position missing in UpdatedSampleList):\n")
  print(missing_join[, c("position","mouseID","treatment","timePoint")])
}

missing_parse <- merged[!is.na(merged$ageSex) & (is.na(merged$age) | is.na(merged$sex)), , drop = FALSE]
if (nrow(missing_parse) > 0) {
  cat("\nRows with ageSex present but unparsable into age/sex:\n")
  print(missing_parse[, c("position","ageSex")])
}

dup_ids <- merged$sampleID[duplicated(merged$sampleID)]
if (length(dup_ids) > 0) {
  cat("\nDuplicate sampleID detected:\n")
  print(unique(dup_ids))
}

cat("\nTotal rows:", nrow(merged), "\n")
cat("===============================================\n\n")

final <- merged[, c("position","mouseID","treatment","timePoint","ageSex","age","sex","sampleID")]
fwrite(final, out_path, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote:\n  ", out_path, "\n", sep = "")
