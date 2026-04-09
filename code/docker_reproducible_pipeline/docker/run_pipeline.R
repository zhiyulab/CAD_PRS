#!/usr/bin/env Rscript
###############################################################################
# Unified Primary Analysis Pipeline
# Combines all 15 R scripts from the primary analysis into one sequential script.
#
# Usage:
#   In Docker:  Rscript run_pipeline.R
#   Locally:    Rscript run_pipeline.R --data_dir=/path/to/data --output_dir=/path/to/output
#
# DATA_DIR should contain:
#   ukbb_pheno/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt
#   intermediate_cad/20250615_pheno_noGP.tsv.gz
#   le8_scores/updated_LE8score_wComposite.txt
#   prs_58/  (62 .txt.gz files, 4 excluded)
#   prs_lipid_removed_3/  (3 .txt files)
#   sliang/pce_zy_2.tsv
#   sliang/ukb_exome_450k_fh.carrier.txt
#   sliang/family_history_final.csv
#   sliang/cohort_p_sex_specific.csv
#   sliang/SMuRFless_CAD_under55.csv
#   sliang/UKB_PRS.csv
###############################################################################

# === Parse command-line arguments ===
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- "/data"
OUTPUT_DIR <- "/output"
RUN_FROM <- 1      # start from this step
RUN_ONLY <- NULL   # if set, only run this step
FRESH <- FALSE     # if TRUE, clear caches

for (arg in args) {
  if (grepl("^--data_dir=", arg))   DATA_DIR   <- sub("^--data_dir=", "", arg)
  if (grepl("^--output_dir=", arg)) OUTPUT_DIR <- sub("^--output_dir=", "", arg)
  if (grepl("^--from=", arg))       RUN_FROM   <- as.numeric(sub("^--from=", "", arg))
  if (grepl("^--only=", arg))       RUN_ONLY   <- as.numeric(sub("^--only=", "", arg))
  if (arg == "--fresh")             FRESH      <- TRUE
}

cat("=== Pipeline Configuration ===\n")
cat("DATA_DIR:  ", DATA_DIR, "\n")
cat("OUTPUT_DIR:", OUTPUT_DIR, "\n")
if (!is.null(RUN_ONLY)) { cat("RUN ONLY:   Step", RUN_ONLY, "\n")
} else if (RUN_FROM > 1) { cat("RUN FROM:   Step", RUN_FROM, "\n") }
cat("\n")

# Step numbering: 1=curation, 2=table one, 3=PRS pickup (3a-3f), 4=regression,
#                 5=comparing, 6=risk factors (6a-6c), 7=plots (7a-7b)
should_run <- function(step_num) {
  if (!is.null(RUN_ONLY)) return(step_num == RUN_ONLY)
  return(step_num >= RUN_FROM)
}

# === Load all required libraries ===
library(data.table)
library(readr)
library(tableone)
library(boot)
library(speedglm)
library(purrr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggh4x)
library(forcats)

# === Set random seed for reproducibility (bootstrap CIs, etc.) ===
set.seed(42)

# === Clear caches if --fresh ===
if (FRESH) {
  cat("Clearing all caches...\n")
  unlink(file.path(OUTPUT_DIR, ".cache"), recursive = TRUE)
}

# === Timing and resource tracking ===
pipeline_start <- proc.time()
step_times <- list()
track_step <- function(name) {
  step_times[[name]] <<- proc.time()
}
report_step <- function(name) {
  elapsed <- (proc.time() - step_times[[name]])["elapsed"]
  mem_mb <- as.numeric(gc(verbose = FALSE)[2, 2])  # max used MB
  cat(sprintf("  [%s] %.1f sec | Peak RAM: %.0f MB\n", name, elapsed, mem_mb))
}

# === Create output subdirectories ===
dir.create(file.path(OUTPUT_DIR, "1_curation"),              recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "2_table_one"),             recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "3_prs_pickup/benchmark"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "3_prs_pickup/SA_lipid_removed"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "3_prs_pickup/SA_benchmark_dist"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "3_prs_pickup/SA_precision_recall"),recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "3_prs_pickup/SA_nonbenchmark"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "4_regression"),            recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "5_comparing_individual"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "6_comparing_risk_factor"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "7_plots"),                 recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Load and QC phenotype data (with caching)
###############################################################################
CACHE_DIR <- file.path(OUTPUT_DIR, ".cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

load_pheno_qc <- function() {
  cache_file <- file.path(CACHE_DIR, "pheno_qc.rds")
  if (file.exists(cache_file)) {
    cat("  Loading phenotype data from cache...\n")
    return(readRDS(cache_file))
  }
  cat("  Loading phenotype data (first run, caching for next time)...\n")
  pheno <- fread(file.path(DATA_DIR, "ukbb_pheno", "ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt"))
  pheno$IDs_toRemove_SampleQCv3 <- ifelse(
    (!(pheno$Submitted_Gender == pheno$Inferred_Gender) | pheno$Non_Consented == 1),
    1, 0
  )
  pheno <- pheno[pheno$IDs_toRemove_SampleQCv3 == 0, ]
  pheno <- as.data.frame(pheno)
  saveRDS(pheno, cache_file)
  return(pheno)
}

###############################################################################
# Read 58 PRS files, build UKB_PRS_58 with cutoffs
###############################################################################
load_prs_58 <- function() {
  cache_file <- file.path(CACHE_DIR, "UKB_PRS_58.rds")
  cache_map <- file.path(CACHE_DIR, "prs_name_map.rds")
  if (file.exists(cache_file) && file.exists(cache_map)) {
    cat("  Loading 58 PRS from cache...\n")
    prs_name_map <<- readRDS(cache_map)
    return(readRDS(cache_file))
  }
  cat("  Loading 58 PRS files (first run, caching for next time)...\n")
  prs_dir <- file.path(DATA_DIR, "prs_58")
  all_files <- list.files(path = prs_dir, pattern = "^PC4\\.(.*)\\.adjNormPRS\\.txt\\.gz$", full.names = TRUE)
  excluded_ids <- c("PGS000116", "PGS004899", "PGS004888", "PGS004237")
  filtered_files <- all_files[!grepl(paste(excluded_ids, collapse = "|"), all_files)]
  filtered_files <- sort(filtered_files)

  # Extract PGS IDs from filenames for labeling
  pgs_ids <- gsub(".*PC4\\.(PGS[0-9]+)\\.adjNormPRS.*", "\\1", basename(filtered_files))
  prs_name_map <<- data.frame(df_name = paste0("df", seq_along(filtered_files)), pgs_id = pgs_ids,
                                stringsAsFactors = FALSE)

  data_list <- lapply(filtered_files, function(f) read.table(f, header = TRUE, sep = "\t"))
  n_prs <- length(data_list)
  cat("    Found", n_prs, "PRS files after exclusions\n")

  # Assign to df1..dfN in calling environment
  for (i in seq_along(data_list)) {
    assign(paste0("df", i), data_list[[i]], envir = parent.frame())
  }

  # Build combined PRS matrix
  prs_list <- lapply(data_list, function(d) d[[4]])
  UKB_PRS_58 <- as.data.frame(prs_list)
  colnames(UKB_PRS_58) <- paste0("df", 1:n_prs)
  UKB_PRS_58$id <- data_list[[1]]$IID
  UKB_PRS_58 <- UKB_PRS_58[, c(ncol(UKB_PRS_58), 1:(ncol(UKB_PRS_58) - 1))]

  # Add top-10% cutoffs
  for (i in 1:n_prs) {
    col_name <- paste0("df", i)
    cutoff_val <- quantile(UKB_PRS_58[[col_name]], 0.9, na.rm = TRUE)
    UKB_PRS_58[[paste0(col_name, "_cutoff")]] <- ifelse(UKB_PRS_58[[col_name]] >= cutoff_val, 1, 0)
  }

  saveRDS(UKB_PRS_58, cache_file)
  saveRDS(prs_name_map, cache_map)
  return(UKB_PRS_58)
}

###############################################################################
# Bootstrap capture proportion for cutoff columns
###############################################################################
boot_capture <- function(cohort_df, n_scores, prefix = "df", R = 1000) {
  results <- lapply(1:n_scores, function(i) {
    col_name <- paste0(prefix, i, "_cutoff")
    x <- cohort_df[[col_name]]
    x <- x[!is.na(x)]

    boot_out <- boot(
      data = x,
      statistic = function(x, indices) mean(x[indices] == 1),
      R = R
    )
    ci <- boot.ci(boot_out, type = "perc")

    data.frame(
      PRS = paste0(prefix, i),
      prop = mean(x == 1),
      low_ci = ci$percent[4],
      high_ci = ci$percent[5]
    )
  })
  do.call(rbind, results)
}

###############################################################################
# Relabel df1..df58 to PGS IDs
###############################################################################
relabel_prs <- function(df, prs_col = "PRS") {
  if (exists("prs_name_map")) {
    m <- prs_name_map
    idx <- match(df[[prs_col]], m$df_name)
    df[[prs_col]] <- ifelse(is.na(idx), df[[prs_col]], m$pgs_id[idx])
  }
  return(df)
}

###############################################################################
# STEP 1: Curation
###############################################################################
# Always load pheno (needed by many steps)
pheno <- load_pheno_qc()

if (should_run(1)) {
cat("\n========== STEP 1: Curation ==========\n")
track_step("Step 1")

cache_intermed <- file.path(CACHE_DIR, "pheno_noGP.rds")
if (file.exists(cache_intermed)) {
  cat("  Loading intermediate CAD pheno from cache...\n")
  X20250615_pheno_noGP <- readRDS(cache_intermed)
} else {
  X20250615_pheno_noGP <- as.data.frame(fread(file.path(DATA_DIR, "intermediate_cad", "20250615_pheno_noGP.tsv.gz")))
  saveRDS(X20250615_pheno_noGP, cache_intermed)
}
pce_zy_2 <- as.data.frame(read_delim(file.path(DATA_DIR, "sliang", "pce_zy_2.tsv"), delim = "\t", escape_double = FALSE, trim_ws = TRUE))
PCE <- pce_zy_2[, c(1, 12)]

# Select individuals with intermediate CAD
ID_intermediate <- subset(X20250615_pheno_noGP, has_disease_Coronary_Artery_Disease_INTERMEDIATE == 1)
subset_data_intermediate <- pheno[pheno$id %in% ID_intermediate$sample_id, ]

# Select key phenotype columns and remove missing
data_intermediate <- subset_data_intermediate[, c(1, 14, 15, 23, 22, 74, 184, 178, 187, 2226, 2227, 2228)]
data_intermediate <- merge(data_intermediate, PCE, by = "id", all.x = TRUE)
data_intermediate <- data_intermediate[complete.cases(data_intermediate), ]
colnames(data_intermediate)[ncol(data_intermediate)] <- "PCE"

# Convert HbA1c to NGSP units
data_intermediate$HbA1c <- data_intermediate$Glycated.haemoglobin..HbA1c. * 0.09148 + 2.152

# Remove individuals with SMuRFs
data_intermediate$IDs_toRemove <- ifelse(
  data_intermediate$HbA1c >= 6.5 |
    data_intermediate$PCE >= 0.05 |
    data_intermediate$LDL.direct >= 4.9 |
    data_intermediate$Diabetes_Med == 1 |
    data_intermediate$LipidLowering_Med == 1,
  1, 0
)
cohort_p <- subset(data_intermediate, IDs_toRemove == 0)

# Calculate onset age and restrict to early onset <55
time <- X20250615_pheno_noGP[X20250615_pheno_noGP$sample_id %in% cohort_p$id, ]
time <- time[, c(1, 178, 179, 180, 181)]
colnames(time)[1] <- "id"

cohort_p_all <- merge(time, cohort_p, by = "id", all.x = TRUE)
cohort_p_all$onset_age <- cohort_p_all$age + cohort_p_all$followup_Coronary_Artery_Disease_INTERMEDIATE
cohort_p_55 <- subset(cohort_p_all, onset_age < 55)

cohort_p_all <- as.data.frame(cohort_p_all)
cohort_p_55 <- as.data.frame(cohort_p_55)
fwrite(cohort_p_all, file.path(OUTPUT_DIR, "1_curation", "cohort_p_all.csv"))
fwrite(cohort_p_55, file.path(OUTPUT_DIR, "1_curation", "cohort_p_55.csv"))

complete_data_intermediate <- subset_data_intermediate[subset_data_intermediate$id %in% data_intermediate$id, ]
not_cohort_p <- as.data.frame(complete_data_intermediate[!complete_data_intermediate$id %in% cohort_p_all$id, ])
fwrite(not_cohort_p, file.path(OUTPUT_DIR, "1_curation", "not_cohort_p.csv"))

cat("  Cohort_p_all:", nrow(cohort_p_all), "individuals\n")
cat("  Cohort_p_55:", nrow(cohort_p_55), "individuals\n")
cat("  Not_cohort_p:", nrow(not_cohort_p), "individuals\n")
# Save key variables for step-skipping
saveRDS(list(cohort_p_all=cohort_p_all, cohort_p_55=cohort_p_55, not_cohort_p=not_cohort_p,
             X20250615_pheno_noGP=X20250615_pheno_noGP, pce_zy_2=pce_zy_2, PCE=PCE),
        file.path(CACHE_DIR, "step1_vars.rds"))
report_step("Step 1")

} else {
  cat("\n========== STEP 1: Skipped (loading cached) ==========\n")
  s1 <- readRDS(file.path(CACHE_DIR, "step1_vars.rds"))
  cohort_p_all <- s1$cohort_p_all; cohort_p_55 <- s1$cohort_p_55; not_cohort_p <- s1$not_cohort_p
  X20250615_pheno_noGP <- s1$X20250615_pheno_noGP; pce_zy_2 <- s1$pce_zy_2; PCE <- s1$PCE
  rm(s1)
}

###############################################################################
# STEP 2: Table One
###############################################################################
if (should_run(2)) {
cat("\n========== STEP 2: Table One ==========\n")
track_step("Step 2")

LE8score_composite_ <- read_csv(file.path(DATA_DIR, "le8_scores", "updated_LE8score_wComposite.txt"))
family <- read_csv(file.path(DATA_DIR, "sliang", "family_history_final.csv"))
fh <- read_delim(file.path(DATA_DIR, "sliang", "ukb_exome_450k_fh.carrier.txt"),
                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)
fh <- fh[, c(1, 3)]
colnames(fh)[1] <- "id"

# Use id-only versions for merging (keep as data.frame)
cohort_p_all_ids <- data.frame(id = cohort_p_all$id)
cohort_p_55_ids <- data.frame(id = cohort_p_55$id)
not_cohort_p_ids <- data.frame(id = not_cohort_p$id)

build_tableone_df <- function(id_df) {
  tbl <- merge(id_df, pheno, by = "id", all.x = TRUE)
  tbl <- tbl[, c(1, 12, 14, 74, 180, 2309)]
  tbl <- merge(tbl, PCE, by = "id", all.x = TRUE)
  tbl <- merge(tbl, LE8score_composite_, by = "id", all.x = TRUE)
  tbl <- merge(tbl, family, by = "id", all.x = TRUE)
  tbl <- merge(tbl, fh, by = "id", all.x = TRUE)
  tbl$Sex.y <- NULL
  tbl$age.y <- NULL
  colnames(tbl) <- c("id", "Sex", "Age", "White British Ancestry", "C Reactive Protein",
                      "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                      "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score",
                      "LE8 Sleep Health Points", "LE8 BMI Points",
                      "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score",
                      "LE8 Composite Score", "Family Heart Disease History",
                      "Familial Hypercholesterolemia Variant Status")
  return(tbl)
}

tableone_all <- build_tableone_df(cohort_p_all_ids)
tableone_55 <- build_tableone_df(cohort_p_55_ids)
tableone_not <- build_tableone_df(not_cohort_p_ids)

myVars <- c("Sex", "Age", "White British Ancestry", "C Reactive Protein",
            "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
            "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score",
            "LE8 Sleep Health Points", "LE8 BMI Points",
            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score",
            "LE8 Composite Score", "Family Heart Disease History",
            "Familial Hypercholesterolemia Variant Status")
catVars <- c("Sex", "White British Ancestry", "Family Heart Disease History",
             "Familial Hypercholesterolemia Variant Status")

tab1 <- CreateTableOne(vars = myVars, data = tableone_all, factorVars = catVars)
table1_pall <- print(tab1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_pall, file = file.path(OUTPUT_DIR, "2_table_one", "Table1_pall.csv"))

tab2 <- CreateTableOne(vars = myVars, data = tableone_55, factorVars = catVars)
table1_p55 <- print(tab2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_p55, file = file.path(OUTPUT_DIR, "2_table_one", "Table1_p55.csv"))

tab3 <- CreateTableOne(vars = myVars, data = tableone_not, factorVars = catVars)
table1_pnot <- print(tab3, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_pnot, file = file.path(OUTPUT_DIR, "2_table_one", "table1_pnot.csv"))

# P-values for continuous variables (all vs not)
continuous_cols <- c("Age", "C Reactive Protein", "Townsend Deprivation Index",
                     "Pooled Cohort Equation Risk", "LE8 Diet Score", "LE8 PA Score",
                     "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                     "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score",
                     "LE8 Composite Score")
p_values <- sapply(continuous_cols, function(col) t.test(tableone_all[[col]], tableone_not[[col]])$p.value)
p_values_formatted <- format.pval(p_values, eps = 1e-300)
pval_continuous_risk <- data.frame(variable = continuous_cols, p_value = p_values_formatted)
write.csv(pval_continuous_risk, file = file.path(OUTPUT_DIR, "2_table_one", "pval_continuous_risk.csv"), row.names = FALSE)

# P-values for continuous variables (55 vs all)
p_values <- sapply(continuous_cols, function(col) t.test(tableone_55[[col]], tableone_all[[col]])$p.value)
pval_continuous_age <- data.frame(variable = continuous_cols, p_value = as.numeric(p_values))
write.csv(pval_continuous_age, file = file.path(OUTPUT_DIR, "2_table_one", "pval_continuous_page.csv"), row.names = FALSE)

# P-values for categorical variables (all vs not)
categorical_cols <- c("Sex", "White British Ancestry", "Family Heart Disease History",
                      "Familial Hypercholesterolemia Variant Status")
get_pvals_risk <- function(col) {
  df_tmp <- rbind(
    data.frame(Group = "highrisk", Var = tableone_all[[col]]),
    data.frame(Group = "nothighrisk", Var = tableone_not[[col]])
  )
  table_data <- table(df_tmp$Group, df_tmp$Var)
  fisher_p <- tryCatch(fisher.test(table_data)$p.value, error = function(e) NA)
  chisq_p <- tryCatch(chisq.test(table_data)$p.value, error = function(e) NA)
  return(c(chi_sq_p = chisq_p, fisher_p = fisher_p))
}
pval_results <- t(sapply(categorical_cols, get_pvals_risk))
pval_categorical_risk <- data.frame(
  variable = rownames(pval_results),
  chi_sq_p = pval_results[, "chi_sq_p"],
  fisher_p = pval_results[, "fisher_p"],
  row.names = NULL
)
write.csv(pval_categorical_risk, file = file.path(OUTPUT_DIR, "2_table_one", "pval_categorical_risk.csv"))

# P-values for categorical variables (all vs 55)
get_pvals_age <- function(col) {
  df_tmp <- rbind(
    data.frame(Group = "all", Var = tableone_all[[col]]),
    data.frame(Group = "earlyonset", Var = tableone_55[[col]])
  )
  table_data <- table(df_tmp$Group, df_tmp$Var)
  fisher_p <- tryCatch(fisher.test(table_data)$p.value, error = function(e) NA)
  chisq_p <- tryCatch(chisq.test(table_data)$p.value, error = function(e) NA)
  return(c(chi_sq_p = chisq_p, fisher_p = fisher_p))
}
pval_results <- t(sapply(categorical_cols, get_pvals_age))
pval_categorical_age <- data.frame(
  variable = rownames(pval_results),
  chi_sq_p = pval_results[, "chi_sq_p"],
  fisher_p = pval_results[, "fisher_p"],
  row.names = NULL
)
write.csv(pval_categorical_age, file = file.path(OUTPUT_DIR, "2_table_one", "pval_categorical_page.csv"))

saveRDS(list(LE8score_composite_=LE8score_composite_, family=family, fh=fh),
        file.path(CACHE_DIR, "step2_vars.rds"))
report_step("Step 2")

} else {
  cat("\n========== STEP 2: Skipped (loading cached) ==========\n")
  s2 <- readRDS(file.path(CACHE_DIR, "step2_vars.rds"))
  LE8score_composite_ <- s2$LE8score_composite_; family <- s2$family; fh <- s2$fh
  rm(s2)
}

###############################################################################
# STEP 3a: PRS Pickup — Race Stratified (benchmark)
###############################################################################
# Always load PRS (needed by steps 3-7)
UKB_PRS_58 <- load_prs_58()

if (should_run(3)) {
cat("\n========== STEP 3a: PRS Pickup (Race Stratified) ==========\n")
track_step("Step 3a")

# Merge with cohort_p_55
cohort_p_race <- cohort_p_55[, c("id", "in_white_British_ancestry_subset")]
cohort_p_race <- merge(cohort_p_race, UKB_PRS_58, by = "id", all.x = TRUE)
cohort_p_other <- cohort_p_race[cohort_p_race$in_white_British_ancestry_subset == 0, ]
cohort_p_white <- cohort_p_race[cohort_p_race$in_white_British_ancestry_subset == 1, ]

# Bootstrap: all
cat("  Bootstrap: all ancestries...\n")
boot_race_all <- boot_capture(cohort_p_race, 58)
fwrite(boot_race_all, file.path(OUTPUT_DIR, "3_prs_pickup/benchmark", "boot_summary_df_55.csv"))

# Bootstrap: white
cat("  Bootstrap: white ancestry...\n")
boot_race_white <- boot_capture(cohort_p_white, 58)
fwrite(boot_race_white, file.path(OUTPUT_DIR, "3_prs_pickup/benchmark", "boot_summary_df_55_white.csv"))

# Bootstrap: other
cat("  Bootstrap: other ancestry...\n")
boot_race_other <- boot_capture(cohort_p_other, 58)
fwrite(boot_race_other, file.path(OUTPUT_DIR, "3_prs_pickup/benchmark", "boot_summary_df_55_other.csv"))

# Code after stop(): calculate pickup sums
cutoff_cols <- grep("_cutoff$", names(cohort_p_race), value = TRUE)
cohort_p_sum <- cohort_p_race[, c("id", cutoff_cols), drop = FALSE]
cohort_p_sum$sum <- rowSums(cohort_p_sum[, cutoff_cols, drop = FALSE], na.rm = TRUE)
fwrite(cohort_p_sum, file.path(OUTPUT_DIR, "1_curation", "cohort_p_sum.csv"))

# Build combined boot_summary_df_55 for plotting (race stratified)
# This is what capture_plot.R expects
boot_race_all$Group <- "Total"
boot_race_white$Group <- "British_White"
boot_race_other$Group <- "other"
boot_summary_df_55 <- rbind(boot_race_all, boot_race_white, boot_race_other)
boot_summary_df_55 <- boot_summary_df_55[, c("PRS", "Group", "prop", "low_ci", "high_ci")]

report_step("Step 3a")

###############################################################################
# STEP 3b: PRS Pickup — Sex Stratified (benchmark)
###############################################################################
cat("\n========== STEP 3b: PRS Pickup (Sex Stratified) ==========\n")
track_step("Step 3b")

cohort_p_sex_specific <- read.csv(file.path(DATA_DIR, "sliang", "cohort_p_sex_specific.csv"))
cohort_p_sex <- cohort_p_sex_specific[, c("id", "Sex")]
cohort_p_sex <- merge(cohort_p_sex, UKB_PRS_58, by = "id", all.x = TRUE)

cohort_p_m <- cohort_p_sex[cohort_p_sex$Sex == "Male", ]
cohort_p_f <- cohort_p_sex[cohort_p_sex$Sex == "Female", ]

fwrite(cohort_p_m, file.path(OUTPUT_DIR, "3_prs_pickup/benchmark", "cohort_p_m.csv"))
fwrite(cohort_p_f, file.path(OUTPUT_DIR, "3_prs_pickup/benchmark", "cohort_p_f.csv"))

cat("  Bootstrap: male...\n")
boot_sex_m <- boot_capture(cohort_p_m, 58)
fwrite(boot_sex_m, file.path(OUTPUT_DIR, "3_prs_pickup/benchmark", "boot_summary_df_55_m.csv"))

cat("  Bootstrap: female...\n")
boot_sex_f <- boot_capture(cohort_p_f, 58)
fwrite(boot_sex_f, file.path(OUTPUT_DIR, "3_prs_pickup/benchmark", "boot_summary_df_55_f.csv"))

# Build combined boot_summary_df_55_sex for plotting
# Re-use boot_race_all as "Total" for sex plot (same overall cohort)
boot_sex_total <- boot_race_all
boot_sex_total$Group <- "Total"
boot_sex_m$Group <- "Male"
boot_sex_f$Group <- "Female"
boot_summary_df_55_sex <- rbind(boot_sex_total, boot_sex_m, boot_sex_f)
boot_summary_df_55_sex <- boot_summary_df_55_sex[, c("PRS", "Group", "prop", "low_ci", "high_ci")]

report_step("Step 3b")

###############################################################################
# STEP 3c: SA — Lipid Removed PRS
###############################################################################
cat("\n========== STEP 3c: Lipid Removed PRS ==========\n")
track_step("Step 3c")

lipid_dir <- file.path(DATA_DIR, "prs_lipid_removed_3")
lipid_files <- list.files(path = lipid_dir, pattern = "^PC4\\.(.*)\\.adjNormPRS\\.txt$", full.names = TRUE)
lipid_files <- sort(lipid_files)
lipid_data_list <- lapply(lipid_files, function(f) read.table(f, header = TRUE, sep = "\t"))

prs_list_lipid <- lapply(lipid_data_list, function(d) d[[4]])
lipid_removed_prs <- as.data.frame(prs_list_lipid)
colnames(lipid_removed_prs) <- paste0("df", 1:3)
lipid_removed_prs$id <- lipid_data_list[[1]]$IID
lipid_removed_prs <- lipid_removed_prs[, c(ncol(lipid_removed_prs), 1:(ncol(lipid_removed_prs) - 1))]

for (i in 1:3) {
  col_name <- paste0("df", i)
  cutoff_val <- quantile(lipid_removed_prs[[col_name]], 0.9, na.rm = TRUE)
  lipid_removed_prs[[paste0(col_name, "_cutoff")]] <- ifelse(lipid_removed_prs[[col_name]] >= cutoff_val, 1, 0)
}
fwrite(lipid_removed_prs, file.path(OUTPUT_DIR, "3_prs_pickup/SA_lipid_removed", "lipid_removed_prs.csv"))

cohort_p_lipid <- cohort_p_55[, c(1, 10)]
cohort_p_lipid <- merge(cohort_p_lipid, lipid_removed_prs, by = "id", all.x = TRUE)
cohort_p_lipid_other <- cohort_p_lipid[cohort_p_lipid$in_white_British_ancestry_subset == 0, ]
cohort_p_lipid_white <- cohort_p_lipid[cohort_p_lipid$in_white_British_ancestry_subset == 1, ]

fwrite(cohort_p_lipid, file.path(OUTPUT_DIR, "3_prs_pickup/SA_lipid_removed", "cohort_p.csv"))
fwrite(cohort_p_lipid_other, file.path(OUTPUT_DIR, "3_prs_pickup/SA_lipid_removed", "cohort_p_other.csv"))
fwrite(cohort_p_lipid_white, file.path(OUTPUT_DIR, "3_prs_pickup/SA_lipid_removed", "cohort_p_white.csv"))

cat("  Bootstrap: lipid removed all...\n")
boot_lipid_all <- boot_capture(cohort_p_lipid, 3)
fwrite(boot_lipid_all, file.path(OUTPUT_DIR, "3_prs_pickup/SA_lipid_removed", "boot_summary_df_55.csv"))

cat("  Bootstrap: lipid removed white...\n")
boot_lipid_white <- boot_capture(cohort_p_lipid_white, 3)
fwrite(boot_lipid_white, file.path(OUTPUT_DIR, "3_prs_pickup/SA_lipid_removed", "boot_summary_df_55_white.csv"))

cat("  Bootstrap: lipid removed other...\n")
boot_lipid_other <- boot_capture(cohort_p_lipid_other, 3)
fwrite(boot_lipid_other, file.path(OUTPUT_DIR, "3_prs_pickup/SA_lipid_removed", "boot_summary_df_55_other.csv"))

report_step("Step 3c")

###############################################################################
# STEP 3d: SA — Benchmark Distribution
###############################################################################
cat("\n========== STEP 3d: Benchmark Distribution ==========\n")
track_step("Step 3d")

# Use UKB_PRS_58 (already in memory) as UKB_PRS equivalent
# The script uses UKB_PRS$df53 and cohort_p_55$id
UKB_PRS_dist <- UKB_PRS_58  # full population PRS

probs <- seq(0, 1, by = 0.1)
common_counts <- numeric(length(probs) - 1)
bin_labels <- character(length(probs) - 1)

for (i in 1:(length(probs) - 1)) {
  lower <- quantile(UKB_PRS_dist$df53, probs[i], na.rm = TRUE)
  upper <- quantile(UKB_PRS_dist$df53, probs[i + 1], na.rm = TRUE)
  bin_ids <- UKB_PRS_dist$id[UKB_PRS_dist$df53 >= lower & UKB_PRS_dist$df53 < upper]
  common_counts[i] <- sum(bin_ids %in% cohort_p_55$id)
  bin_labels[i] <- paste0(probs[i] * 100, "-", probs[i + 1] * 100, "%")
}

df_plot <- data.frame(
  bin = factor(bin_labels, levels = rev(bin_labels)),
  common_count = common_counts
)

p_dist <- ggplot(df_plot, aes(x = bin, y = common_count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "PRS Percentile Bin", y = "Number of Common IDs",
       title = "Number of IDs in cohort by PRS Percentile Bin") +
  theme_minimal() +
  coord_flip()

ggsave(file.path(OUTPUT_DIR, "3_prs_pickup/SA_benchmark_dist", "benchmark_distribution.png"),
       plot = p_dist, width = 6, height = 4, dpi = 300)

report_step("Step 3d")

###############################################################################
# STEP 3e: SA — Precision Recall
###############################################################################
cat("\n========== STEP 3e: Precision Recall ==========\n")
track_step("Step 3e")

UKB_PRS_pr <- UKB_PRS_58[, c("id", "df53")]

cutoff_val <- quantile(UKB_PRS_pr$df53, 0.9, na.rm = TRUE)
UKB_PRS_pr$cutoff_10 <- ifelse(UKB_PRS_pr$df53 >= cutoff_val, 1, 0)
cutoff_val <- quantile(UKB_PRS_pr$df53, 0.8, na.rm = TRUE)
UKB_PRS_pr$cutoff_20 <- ifelse(UKB_PRS_pr$df53 >= cutoff_val, 1, 0)
cutoff_val <- quantile(UKB_PRS_pr$df53, 0.99, na.rm = TRUE)
UKB_PRS_pr$cutoff_1 <- ifelse(UKB_PRS_pr$df53 >= cutoff_val, 1, 0)
cutoff_val <- quantile(UKB_PRS_pr$df53, 0.95, na.rm = TRUE)
UKB_PRS_pr$cutoff_5 <- ifelse(UKB_PRS_pr$df53 >= cutoff_val, 1, 0)

colnames(UKB_PRS_pr) <- c("id", "df53", "10% cutoff", "20% cutoff", "1% cutoff", "5% cutoff")
UKB_PRS_pr$case <- ifelse(UKB_PRS_pr$id %in% cohort_p_55$id, 1, 0)

threshold_cols <- c("1% cutoff", "5% cutoff", "10% cutoff", "20% cutoff")

pr_df <- UKB_PRS_pr %>%
  pivot_longer(cols = all_of(threshold_cols),
               names_to = "threshold",
               values_to = "pred") %>%
  group_by(threshold) %>%
  summarise(
    TP = sum(pred == 1 & case == 1),
    FP = sum(pred == 1 & case == 0),
    FN = sum(pred == 0 & case == 1),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    .groups = "drop"
  )

pr_plot <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_line(color = "grey40") +
  geom_point(size = 3, color = "steelblue") +
  geom_text_repel(
    aes(label = threshold),
    size = 3.5, box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "Recall", y = "Precision")

ggsave(file.path(OUTPUT_DIR, "3_prs_pickup/SA_precision_recall", "precision_recall_curve.png"),
       plot = pr_plot, width = 6, height = 4, dpi = 300)

report_step("Step 3e")

###############################################################################
# STEP 3f: SA — Non-benchmark (not_cohort_p)
###############################################################################
cat("\n========== STEP 3f: Non-benchmark PRS Pickup ==========\n")
track_step("Step 3f")

cohort_p_nonbench <- not_cohort_p[, c("id", "in_white_British_ancestry_subset")]
cohort_p_nonbench <- merge(cohort_p_nonbench, UKB_PRS_58, by = "id", all.x = TRUE)
cohort_p_nonbench_other <- cohort_p_nonbench[cohort_p_nonbench$in_white_British_ancestry_subset == 0, ]
cohort_p_nonbench_white <- cohort_p_nonbench[cohort_p_nonbench$in_white_British_ancestry_subset == 1, ]

cat("  Bootstrap: nonbenchmark all...\n")
boot_nonbench_all <- boot_capture(cohort_p_nonbench, 58)
fwrite(boot_nonbench_all, file.path(OUTPUT_DIR, "3_prs_pickup/SA_nonbenchmark", "boot_summary_df_55.csv"))

cat("  Bootstrap: nonbenchmark white...\n")
boot_nonbench_white <- boot_capture(cohort_p_nonbench_white, 58)
fwrite(boot_nonbench_white, file.path(OUTPUT_DIR, "3_prs_pickup/SA_nonbenchmark", "boot_summary_df_55_white.csv"))

cat("  Bootstrap: nonbenchmark other...\n")
boot_nonbench_other <- boot_capture(cohort_p_nonbench_other, 58)
fwrite(boot_nonbench_other, file.path(OUTPUT_DIR, "3_prs_pickup/SA_nonbenchmark", "boot_summary_df_55_other.csv"))

saveRDS(list(boot_race_all=boot_race_all, boot_summary_df_55=boot_summary_df_55,
             boot_summary_df_55_sex=boot_summary_df_55_sex, cohort_p_race=cohort_p_race),
        file.path(CACHE_DIR, "step3_vars.rds"))
report_step("Step 3f")

} else {
  cat("\n========== STEP 3: Skipped (loading cached) ==========\n")
  s3 <- readRDS(file.path(CACHE_DIR, "step3_vars.rds"))
  boot_race_all <- s3$boot_race_all; boot_summary_df_55 <- s3$boot_summary_df_55
  boot_summary_df_55_sex <- s3$boot_summary_df_55_sex; cohort_p_race <- s3$cohort_p_race
  rm(s3)
}

###############################################################################
# STEP 4: Regression (OR)
###############################################################################
if (should_run(4)) {
cat("\n========== STEP 4: Regression (OR) ==========\n")
track_step("Step 4")

# UKB_PRS.csv is essentially UKB_PRS_58 but re-read from generated data
# Use UKB_PRS_58 directly (id + df1..df58 columns, no cutoff cols)
UKB_PRS_for_reg <- UKB_PRS_58[, c("id", paste0("df", 1:58))]

pheno1 <- pheno[, c(1:14)]
names(UKB_PRS_for_reg)[1] <- "id"
regression <- merge(data.frame(pheno1), UKB_PRS_for_reg, by = "id", all.y = TRUE)
race_col <- pheno[, c("id", "in_white_British_ancestry_subset")]
regression <- merge(regression, race_col, by = "id", all.x = TRUE)
outcome <- X20250615_pheno_noGP[, c(1, 178, 179, 180, 181)]
names(outcome)[1] <- "id"
regression <- merge(regression, data.frame(outcome), by = "id", all.x = TRUE)
write.table(regression, file.path(OUTPUT_DIR, "4_regression", "regression.txt"), sep = "\t", row.names = FALSE)

regression_other <- regression[regression$in_white_British_ancestry_subset == 0, ]
regression_white <- regression[regression$in_white_British_ancestry_subset == 1, ]

exposures <- paste0("df", 1:58)

run_logistic <- function(var, data) {
  formula <- as.formula(paste("has_disease_Coronary_Artery_Disease_INTERMEDIATE ~", var, "+ age + Sex + PC1 + PC2 + PC3 + PC4"))
  model <- speedglm(formula, data = data, family = binomial())
  coefs <- summary(model)$coefficients
  or <- exp(coefs[var, "Estimate"])
  se <- coefs[var, "Std. Error"]
  ci_low <- exp(coefs[var, "Estimate"] - 1.96 * se)
  ci_high <- exp(coefs[var, "Estimate"] + 1.96 * se)
  p <- coefs[var, "Pr(>|z|)"]
  tibble(Exposure = var, OR = or, CI_lower = ci_low, CI_upper = ci_high, P_value = p)
}

cat("  Running logistic regression: all...\n")
results_or_all <- map_dfr(exposures, run_logistic, data = regression)
cat("  Running logistic regression: white...\n")
results_or_white <- map_dfr(exposures, run_logistic, data = regression_white)
cat("  Running logistic regression: other...\n")
results_or_other <- map_dfr(exposures, run_logistic, data = regression_other)

fwrite(results_or_all, file.path(OUTPUT_DIR, "4_regression", "results.csv"))
fwrite(results_or_other, file.path(OUTPUT_DIR, "4_regression", "results_other.csv"))
fwrite(results_or_white, file.path(OUTPUT_DIR, "4_regression", "results_white.csv"))

# Build combined OR results for plotting
results_or_all$Group <- "Total"
results_or_white$Group <- "British_White"
results_or_other$Group <- "other"
results_or_combined <- rbind(results_or_all, results_or_white, results_or_other)
results_or_combined <- results_or_combined[, c("Exposure", "Group", "OR", "CI_lower", "CI_upper")]
colnames(results_or_combined) <- c("PRS", "Group", "OR", "low_ci", "high_ci")

saveRDS(list(results_or_combined=results_or_combined), file.path(CACHE_DIR, "step4_vars.rds"))
report_step("Step 4")

} else {
  cat("\n========== STEP 4: Skipped (loading cached) ==========\n")
  s4 <- readRDS(file.path(CACHE_DIR, "step4_vars.rds"))
  results_or_combined <- s4$results_or_combined
  rm(s4)
}

###############################################################################
# STEP 5: Comparing Individuals by Capture
###############################################################################
if (should_run(5)) {
cat("\n========== STEP 5: Comparing Individuals by Capture ==========\n")
track_step("Step 5")

cohort_p_sum_loaded <- as.data.frame(fread(file.path(OUTPUT_DIR, "1_curation", "cohort_p_sum.csv")))
cohort_p_sum_loaded <- cohort_p_sum_loaded[!is.na(cohort_p_sum_loaded$df1_cutoff), ]
cohort_p_sum_loaded <- cohort_p_sum_loaded[, c("id", "sum")]
onsetage <- data.frame(id = cohort_p_55$id, onset_age = cohort_p_55$onset_age)

ID <- merge(cohort_p_sum_loaded, pheno, by = "id", all.x = TRUE)
ID <- merge(ID, LE8score_composite_, by = "id", all.x = TRUE)
ID <- merge(ID, pce_zy_2, by = "id", all.x = TRUE)
# Use original column indices — these are positions in the merged ID dataframe
# c(1=id, 2=sum, 13=Sex, 15=SmokingStatus, 16=Townsend, 2310+=LE8 scores, 1392=Race, etc.)
# Fall back to name-based selection if indices are out of range
if (ncol(ID) >= 2334) {
  peopleprofile <- ID[, c(1,2,13,15,16,2310,2315,2316,2317,2318,2319,2320,2321,2322,2323,1392,21,647,939,2334,174,189,181)]
} else {
  # Name-based fallback for smaller datasets
  if ("Sex.x" %in% names(ID) && !"Sex" %in% names(ID)) names(ID)[names(ID) == "Sex.x"] <- "Sex"
  if ("age.x" %in% names(ID) && !"age" %in% names(ID)) names(ID)[names(ID) == "age.x"] <- "age"
  if ("SmokingStatus.x" %in% names(ID) && !"SmokingStatus" %in% names(ID)) names(ID)[names(ID) == "SmokingStatus.x"] <- "SmokingStatus"
  pp_col_names <- c("id", "sum", "Sex", "age", "SmokingStatus",
                     "Townsend",
                     "diet_score", "PA_score", "final_smoking_score",
                     "sleep_health_points", "bmi_points",
                     "Blood_lipid_points", "HbA1c_score", "BP_Score",
                     "LE8_composite",
                     "in_white_British_ancestry_subset",
                     "Alcohol_intake_freq",
                     "Chronic_kidney_disease", "Rheumatoid_arthritis",
                     "pce_goff",
                     "Apolipoprotein.B.x", "Lipoprotein.A",
                     "C.reactive.protein")
  pp_col_names <- pp_col_names[pp_col_names %in% names(ID)]
  peopleprofile <- ID[, pp_col_names]
}
cat("  peopleprofile:", nrow(peopleprofile), "rows x", ncol(peopleprofile), "cols\n")
peopleprofile <- merge(peopleprofile, onsetage, by = "id", all.x = TRUE)
peopleprofile <- merge(peopleprofile, family, by = "id", all.x = TRUE)
peopleprofile <- merge(peopleprofile, fh, by = "id", all.x = TRUE)

peopleprofile <- peopleprofile[order(-peopleprofile$sum), ]
rownames(peopleprofile) <- 1:nrow(peopleprofile)

not_picked_up <- peopleprofile[peopleprofile$sum == 0, ]
picked_up <- peopleprofile[peopleprofile$sum > 20, ]

# Fallback for small datasets where thresholds may not apply
if (nrow(picked_up) < 2 || nrow(not_picked_up) < 2) {
  cat("  Adaptive thresholds for small dataset...\n")
  q25 <- quantile(peopleprofile$sum, 0.25, na.rm = TRUE)
  q75 <- quantile(peopleprofile$sum, 0.75, na.rm = TRUE)
  not_picked_up <- peopleprofile[peopleprofile$sum <= q25, ]
  picked_up <- peopleprofile[peopleprofile$sum >= q75, ]
  cat("  Using bottom 25% (sum <=", q25, ") vs top 25% (sum >=", q75, ")\n")
}
if (nrow(not_picked_up) > 0) rownames(not_picked_up) <- 1:nrow(not_picked_up)
if (nrow(picked_up) > 0) rownames(picked_up) <- 1:nrow(picked_up)

# Rename columns for display — map from internal names to labels
col_rename_map <- c(
  "id" = "id", "sum" = "sum",
  "Sex" = "Sex", "Sex.x" = "Sex", "age" = "Age", "age.x" = "Age",
  "SmokingStatus" = "Smoking", "SmokingStatusv2" = "Smoking",
  "Townsend" = "Townsend Deprivation Index", "townsend" = "Townsend Deprivation Index",
  "Race" = "Race",
  "diet_score" = "LE8 Diet Score", "PA_score" = "LE8 PA Score",
  "final_smoking_score" = "LE8 Smoking Score", "sleep_health_points" = "LE8 Sleep Health Points",
  "bmi_points" = "LE8 BMI Points", "Blood_lipid_points" = "LE8 Blood Lipid Points",
  "HbA1c_score" = "LE8 HbA1c Score", "BP_Score" = "LE8 BP Score",
  "LE8_composite" = "LE8 Composite Score",
  "in_white_British_ancestry_subset" = "Race",
  "Alcohol_intake_freq" = "Alcohol intake frequency",
  "Chronic_kidney_disease" = "Chronic Kidney Disease",
  "Rheumatoid_arthritis" = "Rheumatoid arthritis",
  "pce_goff" = "Pooled Cohort Equation Risk",
  "Apolipoprotein.B.x" = "Apolipoprotein B", "Lipoprotein.A" = "Lipoprotein A",
  "C.reactive.protein" = "C Reactive Protein",
  "onset_age" = "Age of CAD Onset",
  "family_history" = "Family Heart Disease History",
  "genotype" = "Familial Hypercholesterolemia Variant Status"
)
# Apply rename for columns that exist
for (old_name in names(col_rename_map)) {
  if (old_name %in% colnames(picked_up)) {
    colnames(picked_up)[colnames(picked_up) == old_name] <- col_rename_map[old_name]
    colnames(not_picked_up)[colnames(not_picked_up) == old_name] <- col_rename_map[old_name]
  }
}

all_desired_vars <- c("Sex", "Age", "Smoking", "Townsend Deprivation Index",
             "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score",
             "LE8 Sleep Health Points", "LE8 BMI Points",
             "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score",
             "LE8 Composite Score", "Race", "Alcohol intake frequency",
             "Chronic Kidney Disease", "Rheumatoid arthritis",
             "Pooled Cohort Equation Risk", "Apolipoprotein B",
             "Lipoprotein A", "Age of CAD Onset", "C Reactive Protein",
             "Family Heart Disease History", "Familial Hypercholesterolemia Variant Status")
myVars5 <- all_desired_vars[all_desired_vars %in% colnames(picked_up)]
all_cat <- c("Sex", "Smoking", "Race", "Alcohol intake frequency",
              "Chronic Kidney Disease", "Rheumatoid arthritis",
              "Family Heart Disease History", "Familial Hypercholesterolemia Variant Status")
catVars5 <- all_cat[all_cat %in% colnames(picked_up)]

if (nrow(not_picked_up) < 2 || nrow(picked_up) < 2) {
  cat("  WARNING: Too few individuals in picked_up or not_picked_up groups. Skipping Step 5 tables/tests.\n")
} else {

tab1_5 <- CreateTableOne(vars = myVars5, data = not_picked_up, factorVars = catVars5)
table1_notpickedup <- print(tab1_5, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE, pDigits = 3)
write.csv(table1_notpickedup, file = file.path(OUTPUT_DIR, "5_comparing_individual", "Table1_notpickedup_p.csv"))

tab2_5 <- CreateTableOne(vars = myVars5, data = picked_up, factorVars = catVars5)
table1_pickedup <- print(tab2_5, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE, pDigits = 3)
write.csv(table1_pickedup, file = file.path(OUTPUT_DIR, "5_comparing_individual", "Table1_pickedup_p.csv"))

all_cont5 <- c("Age", "Townsend Deprivation Index", "LE8 Diet Score",
                      "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points",
                      "LE8 BMI Points", "LE8 Blood Lipid Points", "LE8 HbA1c Score",
                      "LE8 BP Score", "LE8 Composite Score", "Pooled Cohort Equation Risk",
                      "Apolipoprotein B", "Lipoprotein A", "C Reactive Protein", "Age of CAD Onset")
continuous_cols5 <- all_cont5[all_cont5 %in% colnames(picked_up)]
p_values5 <- sapply(continuous_cols5, function(col) tryCatch(t.test(not_picked_up[[col]], picked_up[[col]])$p.value, error = function(e) NA))
pval_continuous5 <- data.frame(variable = continuous_cols5, p_value = as.numeric(p_values5))
write.csv(pval_continuous5, file = file.path(OUTPUT_DIR, "5_comparing_individual", "pval_continuous_p.csv"), row.names = FALSE)

categorical_cols5 <- catVars5
get_pvals5 <- function(col) {
  df_tmp <- rbind(
    data.frame(Group = "Not_Picked_Up", Var = not_picked_up[[col]]),
    data.frame(Group = "Picked_Up", Var = picked_up[[col]])
  )
  table_data <- table(df_tmp$Group, df_tmp$Var)
  fisher_p <- tryCatch(fisher.test(table_data)$p.value, error = function(e) NA)
  chisq_p <- tryCatch(chisq.test(table_data)$p.value, error = function(e) NA)
  return(c(chi_sq_p = chisq_p, fisher_p = fisher_p))
}
pval_results5 <- t(sapply(categorical_cols5, get_pvals5))
pval_categorical5 <- data.frame(
  variable = rownames(pval_results5),
  chi_sq_p = pval_results5[, "chi_sq_p"],
  fisher_p = pval_results5[, "fisher_p"],
  row.names = NULL
)
write.csv(pval_categorical5, file = file.path(OUTPUT_DIR, "5_comparing_individual", "pval_categorical_p.csv"))

} # end if enough individuals

report_step("Step 5")
} else {
  cat("\n========== STEP 5: Skipped ==========\n")
}

###############################################################################
# STEP 6a: Single Metrics (BMI, Lp(a), CRP)
###############################################################################
if (should_run(6)) {
cat("\n========== STEP 6a: Single Metrics ==========\n")
track_step("Step 6a")

# BMI
bmi <- data.frame(pheno[, c(1, 18)])
colnames(bmi) <- c("id", "BMI")
cohort_p_55_ids_only <- data.frame(id = cohort_p_55$id)
cutoff_val <- quantile(bmi$BMI, 0.9, na.rm = TRUE)
bmi$cutoff_10 <- ifelse(bmi$BMI >= cutoff_val, 1, 0)
cohort_p_55_BMI <- merge(cohort_p_55_ids_only, bmi, by = "id", all.x = TRUE)

x <- cohort_p_55_BMI$cutoff_10
x <- x[!is.na(x)]
boot_out <- boot(data = x, statistic = function(x, indices) mean(x[indices] == 1), R = 1000)
ci <- boot.ci(boot_out, type = "perc")
results_bmi <- data.frame(variable = "bmi_cutoff", prop = mean(x == 1),
                           lower_95 = ci$percent[4], upper_95 = ci$percent[5])
fwrite(results_bmi, file.path(OUTPUT_DIR, "6_comparing_risk_factor", "results_bmi.csv"))

# Lp(a) and CRP
lpa_crp <- data.frame(pheno[, c(1, 180, 188)])
colnames(lpa_crp) <- c("id", "crp", "lpa")
lpa_cutoff_val <- quantile(lpa_crp$lpa, 0.9, na.rm = TRUE)
crp_cutoff_val <- quantile(lpa_crp$crp, 0.9, na.rm = TRUE)
lpa_crp$lpa_cutoff <- ifelse(lpa_crp$lpa >= lpa_cutoff_val, 1, 0)
lpa_crp$crp_cutoff <- ifelse(lpa_crp$crp >= crp_cutoff_val, 1, 0)
cohort_p_55_lpa_crp <- merge(cohort_p_55_ids_only, lpa_crp, by = "id", all.x = TRUE)

# Lp(a) bootstrap
x <- cohort_p_55_lpa_crp$lpa_cutoff
x <- x[!is.na(x)]
boot_out <- boot(data = x, statistic = function(x, indices) mean(x[indices] == 1), R = 1000)
ci <- boot.ci(boot_out, type = "perc")
results_lpa <- data.frame(variable = "lpa_cutoff", prop = mean(x == 1),
                           lower_95 = ci$percent[4], upper_95 = ci$percent[5])
fwrite(results_lpa, file.path(OUTPUT_DIR, "6_comparing_risk_factor", "results_lpa.csv"))

# CRP bootstrap
x <- cohort_p_55_lpa_crp$crp_cutoff
x <- x[!is.na(x)]
boot_out <- boot(data = x, statistic = function(x, indices) mean(x[indices] == 1), R = 1000)
ci <- boot.ci(boot_out, type = "perc")
results_crp <- data.frame(variable = "crp_cutoff", prop = mean(x == 1),
                           lower_95 = ci$percent[4], upper_95 = ci$percent[5])
fwrite(results_crp, file.path(OUTPUT_DIR, "6_comparing_risk_factor", "results_crp.csv"))

# Bar plot is generated after step 6b (needs PRS+Lp(a) and PRS+PCE results)

report_step("Step 6a")

###############################################################################
# STEP 6b: Combined Metrics (PRS + Lp(a), PRS + PCE)
###############################################################################
cat("\n========== STEP 6b: Combined Metrics ==========\n")
track_step("Step 6b")

lpa_col <- data.frame(pheno[, c(1, 188)])
colnames(lpa_col) <- c("id", "lpa")
lpa_col <- lpa_col[complete.cases(lpa_col), ]

PCE_comb <- pce_zy_2[, c(1, 12)]
colnames(PCE_comb) <- c("id", "pce")
PCE_comb <- PCE_comb[complete.cases(PCE_comb), ]

UKB_PRS_comb <- UKB_PRS_58[, c("id", paste0("df", 1:58))]
UKB_PRS_comb <- UKB_PRS_comb[complete.cases(UKB_PRS_comb), ]
UKB_PRS_lpa <- merge(UKB_PRS_comb, lpa_col, by = "id")
UKB_PRS_lpa_pce <- merge(UKB_PRS_lpa, PCE_comb, by = "id")

UKB_PRS_lpa_pce$lpa_standardized <- as.numeric(scale(UKB_PRS_lpa_pce$lpa))

for (i in 1:58) {
  col_name <- paste0("df", i)
  UKB_PRS_lpa_pce[[paste0(col_name, "_lpa")]] <- UKB_PRS_lpa_pce[[col_name]] + UKB_PRS_lpa_pce$lpa_standardized
  UKB_PRS_lpa_pce[[paste0(col_name, "_pce")]] <- UKB_PRS_lpa_pce[[col_name]] + UKB_PRS_lpa_pce$pce
}

for (i in 1:58) {
  col_lpa <- paste0("df", i, "_lpa")
  cutoff_val <- quantile(UKB_PRS_lpa_pce[[col_lpa]], 0.9, na.rm = TRUE)
  UKB_PRS_lpa_pce[[paste0(col_lpa, "_cutoff")]] <- ifelse(UKB_PRS_lpa_pce[[col_lpa]] >= cutoff_val, 1, 0)

  col_pce <- paste0("df", i, "_pce")
  cutoff_val <- quantile(UKB_PRS_lpa_pce[[col_pce]], 0.9, na.rm = TRUE)
  UKB_PRS_lpa_pce[[paste0(col_pce, "_cutoff")]] <- ifelse(UKB_PRS_lpa_pce[[col_pce]] >= cutoff_val, 1, 0)
}

cohort_p_55_for_comb <- read.csv(file.path(OUTPUT_DIR, "1_curation", "cohort_p_55.csv"))
UKB_PRS_lpa_pce_p <- UKB_PRS_lpa_pce[UKB_PRS_lpa_pce$id %in% cohort_p_55_for_comb$id, ]
fwrite(UKB_PRS_lpa_pce_p, file.path(OUTPUT_DIR, "6_comparing_risk_factor", "UKB_PRS_lpa_pce_p.csv"))

SMuRFless_55 <- read.csv(file.path(DATA_DIR, "sliang", "SMuRFless_CAD_under55.csv"))
UKB_PRS_lpa_pce_s <- UKB_PRS_lpa_pce[UKB_PRS_lpa_pce$id %in% SMuRFless_55$id, ]
fwrite(UKB_PRS_lpa_pce_s, file.path(OUTPUT_DIR, "6_comparing_risk_factor", "UKB_PRS_lpa_pce_s.csv"))

# Bootstrap PRS+Lpa
x <- UKB_PRS_lpa_pce_p$df53_lpa_cutoff
x <- x[!is.na(x)]
boot_out <- boot(data = x, statistic = function(x, indices) mean(x[indices] == 1), R = 1000)
ci <- boot.ci(boot_out, type = "perc")
results_prs_lpa <- data.frame(variable = "prs_lpa_cutoff", prop = mean(x == 1),
                               lower_95 = ci$percent[4], upper_95 = ci$percent[5])
fwrite(results_prs_lpa, file.path(OUTPUT_DIR, "6_comparing_risk_factor", "results_prs_lpa.csv"))

# Bootstrap PRS+PCE
x <- UKB_PRS_lpa_pce_p$df53_pce_cutoff
x <- x[!is.na(x)]
boot_out <- boot(data = x, statistic = function(x, indices) mean(x[indices] == 1), R = 1000)
ci <- boot.ci(boot_out, type = "perc")
results_prs_pce <- data.frame(variable = "prs_pce_cutoff", prop = mean(x == 1),
                               lower_95 = ci$percent[4], upper_95 = ci$percent[5])
fwrite(results_prs_pce, file.path(OUTPUT_DIR, "6_comparing_risk_factor", "results_prs_pce.csv"))

report_step("Step 6b")

# --- Combined barplot (BMI, CRP, Lp(a), PRS, PRS+Lp(a), PRS+PCE) ---
# Get PRS-only (df53) from step 3a bootstrap
prs_df53 <- boot_race_all[boot_race_all$PRS == "df53", ]
results_prs <- data.frame(variable = "PRS", prop = prs_df53$prop,
                           lower_95 = prs_df53$low_ci, upper_95 = prs_df53$high_ci)

# Combine all metrics and rename for display
all_results <- rbind(results_bmi, results_crp, results_lpa, results_prs,
                      results_prs_lpa, results_prs_pce)
all_results$variable <- c("BMI", "CRP", "Lp(a)", "PRS", "PRS + Lp(a)", "PRS + PCE")
all_results$variable <- factor(all_results$variable,
                                levels = c("BMI", "CRP", "Lp(a)", "PRS", "PRS + Lp(a)", "PRS + PCE"))

p_combined_bar <- ggplot(all_results) +
  geom_bar(aes(x = variable, y = 1), stat = "identity", fill = "grey90") +
  geom_bar(aes(x = variable, y = prop), stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(x = variable, ymin = lower_95, ymax = upper_95), width = 0.2) +
  geom_text(aes(x = variable, y = prop, label = round(prop, 2)), hjust = -0.5, size = 4) +
  coord_flip(clip = "off") +
  theme_minimal(base_size = 12) +
  labs(y = "Proportion (%)", x = "Variable") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(OUTPUT_DIR, "6_comparing_risk_factor", "barplot_with_CI.png"),
       plot = p_combined_bar, width = 6, height = 4, units = "in", dpi = 300, bg = "white")

###############################################################################
# STEP 6c: McNemar Test
###############################################################################
cat("\n========== STEP 6c: McNemar Test ==========\n")
track_step("Step 6c")

# Build the mcnemar dataframe from UKB_PRS_lpa_pce_p
# It needs df53_pce_cutoff and df53_best (which is df53_cutoff from PRS pickup)
mcnemar_df <- UKB_PRS_lpa_pce_p[, c("id", "df53_pce_cutoff")]
# df53_best = top 10% cutoff for df53 from original PRS
prs_cutoff_df53 <- UKB_PRS_58[, c("id", "df53_cutoff")]
colnames(prs_cutoff_df53) <- c("id", "df53_best")
mcnemar_df <- merge(mcnemar_df, prs_cutoff_df53, by = "id", all.x = TRUE)

var1 <- mcnemar_df$df53_pce_cutoff
var2 <- mcnemar_df$df53_best
keep <- complete.cases(var1, var2)
var1_clean <- var1[keep]
var2_clean <- var2[keep]

mcn_result <- mcnemar.test(var1_clean, var2_clean)
cat("  McNemar test result:\n")
print(mcn_result)

# Save result
mcnemar_output <- data.frame(
  statistic = mcn_result$statistic,
  p_value = mcn_result$p.value,
  method = mcn_result$method
)
write.csv(mcnemar_output, file.path(OUTPUT_DIR, "6_comparing_risk_factor", "mcnemar_result.csv"), row.names = FALSE)

report_step("Step 6c")
} else {
  cat("\n========== STEP 6: Skipped ==========\n")
}

###############################################################################
# STEP 7a: Capture Plot (Race and Sex stratified, plus OR)
###############################################################################
if (should_run(7)) {
cat("\n========== STEP 7a: Capture Plot ==========\n")
track_step("Step 7a")

# --- Race-stratified capture plot ---
summary_df_55_final <- relabel_prs(boot_summary_df_55)
colnames(summary_df_55_final) <- c("PRS", "Group", "prop", "low_ci", "high_ci")

summary_df_55_final <- summary_df_55_final %>%
  mutate(
    Group = factor(Group, levels = c("Total", "British_White", "other")),
    Group = case_when(
      Group == "Total" ~ "All Ancestries",
      Group == "British_White" ~ "British White Ancestry",
      Group == "other" ~ "Non-British White Ancestry",
      TRUE ~ as.character(Group)
    )
  ) %>%
  mutate(across(c(prop, low_ci, high_ci), ~ .x * 100))

total_order <- summary_df_55_final %>%
  filter(Group == "All Ancestries") %>%
  arrange(desc(prop)) %>%
  pull(PRS)

summary_df_55_final <- summary_df_55_final %>%
  mutate(PRS = factor(PRS, levels = total_order))

p_race <- ggplot(summary_df_55_final, aes(x = PRS, y = prop, color = Group)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.6,
                position = position_dodge(width = 0.6)) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c(
    "British White Ancestry" = "steelblue",
    "Non-British White Ancestry" = "orange",
    "All Ancestries" = "red"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(y = "Percentage (%)", x = "")

ggsave(file.path(OUTPUT_DIR, "7_plots", "my_plot.png"),
       plot = p_race, width = 12, height = 8, units = "in", dpi = 300)

# --- Sex-stratified capture plot ---
summary_df_55_sex_final <- relabel_prs(boot_summary_df_55_sex)
colnames(summary_df_55_sex_final) <- c("PRS", "Group", "prop", "low_ci", "high_ci")

summary_df_55_sex_final <- summary_df_55_sex_final %>%
  mutate(Group = factor(Group, levels = c("Total", "Male", "Female"))) %>%
  mutate(across(c(prop, low_ci, high_ci), ~ .x * 100))

total_order_sex <- summary_df_55_sex_final %>%
  filter(Group == "Total") %>%
  arrange(desc(prop)) %>%
  pull(PRS)

summary_df_55_sex_final <- summary_df_55_sex_final %>%
  mutate(PRS = factor(PRS, levels = total_order_sex))

p_sex <- ggplot(summary_df_55_sex_final, aes(x = PRS, y = prop, color = Group)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.6,
                position = position_dodge(width = 0.6)) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c(
    "Female" = "steelblue",
    "Male" = "orange",
    "Total" = "red"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(y = "Percentage (%)", x = "")

ggsave(file.path(OUTPUT_DIR, "7_plots", "my_plot_sex.png"),
       plot = p_sex, width = 12, height = 8, units = "in", dpi = 300)

# --- OR plot ---
summary_or_final <- relabel_prs(results_or_combined, prs_col = "PRS")
colnames(summary_or_final) <- c("PRS", "Group", "OR", "low_ci", "high_ci")

summary_or_final <- summary_or_final %>%
  mutate(
    Group = factor(Group, levels = c("Total", "British_White", "other")),
    Group = case_when(
      Group == "Total" ~ "All Ancestries",
      Group == "British_White" ~ "British White Ancestry",
      Group == "other" ~ "Non-British White Ancestry",
      TRUE ~ as.character(Group)
    )
  )

total_order_or <- summary_or_final %>%
  filter(Group == "All Ancestries") %>%
  arrange(desc(OR)) %>%
  pull(PRS)

summary_or_final <- summary_or_final %>%
  mutate(PRS = factor(PRS, levels = total_order_or))

p_or <- ggplot(summary_or_final, aes(x = PRS, y = OR, color = Group)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.6,
                position = position_dodge(width = 0.6)) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c(
    "British White Ancestry" = "steelblue",
    "Non-British White Ancestry" = "orange",
    "All Ancestries" = "red"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(y = "OR/SD", x = "")

ggsave(file.path(OUTPUT_DIR, "7_plots", "my_plot_OR.png"),
       plot = p_or, width = 12, height = 8, units = "in", dpi = 300)

# --- Patchwork: Race-stratified stacked ---
plot_one_group <- function(data, group_name, color,
                           show_x = FALSE, show_y_label = FALSE) {
  ggplot(data, aes(x = PRS, y = prop)) +
    geom_point(size = 3, color = color) +
    geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.6, color = color) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = if (show_x) element_text(angle = 90, hjust = 1) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(),
      axis.ticks.y = element_line(),
      axis.title.y = if (show_y_label) element_text() else element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(title = group_name, y = if (show_y_label) "Percentage (%)" else NULL)
}

p_total <- summary_df_55_final %>%
  filter(Group == "All Ancestries") %>%
  plot_one_group("All Ancestries", "red", show_x = FALSE, show_y_label = FALSE)

p_white <- summary_df_55_final %>%
  filter(Group == "British White Ancestry") %>%
  plot_one_group("British White Ancestry", "steelblue", show_x = FALSE, show_y_label = TRUE)

p_other <- summary_df_55_final %>%
  filter(Group == "Non-British White Ancestry") %>%
  plot_one_group("Non-British White Ancestry", "orange", show_x = TRUE, show_y_label = FALSE)

final_plot_race <- p_total / p_white / p_other +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", color = NA)))

ggsave(file.path(OUTPUT_DIR, "7_plots", "plot_new.png"),
       plot = final_plot_race, width = 12, height = 8, units = "in", dpi = 300)

# --- Patchwork: Sex-stratified stacked ---
p_total_sex <- summary_df_55_sex_final %>%
  filter(Group == "Total") %>%
  plot_one_group("Total", "red", show_x = FALSE, show_y_label = FALSE)

p_female <- summary_df_55_sex_final %>%
  filter(Group == "Female") %>%
  plot_one_group("Female", "steelblue", show_x = FALSE, show_y_label = TRUE)

p_male <- summary_df_55_sex_final %>%
  filter(Group == "Male") %>%
  plot_one_group("Male", "orange", show_x = TRUE, show_y_label = FALSE)

final_plot_sex <- p_total_sex / p_female / p_male +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", color = NA)))

ggsave(file.path(OUTPUT_DIR, "7_plots", "plot_sex_new.png"),
       plot = final_plot_sex, width = 12, height = 8, units = "in", dpi = 300)

report_step("Step 7a")

###############################################################################
# STEP 7b: Slope Plot
###############################################################################
cat("\n========== STEP 7b: Slope Plot ==========\n")
track_step("Step 7b")

# Load the precomputed ranking summary from xlsx (has MGB validation rankings)
library(openxlsx)
prs_ranking_xlsx <- file.path(DATA_DIR, "sliang", "prs_ranking_summary.xlsx")
if (file.exists(prs_ranking_xlsx)) {
  prs_ranking_summary <- read.xlsx(prs_ranking_xlsx)
  ranking <- prs_ranking_summary[, c(2, 9, 10)]
  colnames(ranking) <- c("PRS", "Population", "Individual")
} else {
  # Fallback: compute from pipeline data
  cat("  WARNING: prs_ranking_summary.xlsx not found, computing from pipeline data\n")
  pop_ranking <- relabel_prs(boot_race_all)
  pop_ranking <- pop_ranking[order(-pop_ranking$prop), ]
  pop_ranking$pop_rank <- 1:nrow(pop_ranking)
  individual_props <- sapply(1:58, function(i) {
    col_name <- paste0("df", i, "_cutoff")
    mean(cohort_p_race[[col_name]], na.rm = TRUE)
  })
  ind_ranking <- data.frame(PRS = paste0("df", 1:58), ind_prop = individual_props)
  ind_ranking <- relabel_prs(ind_ranking)
  ind_ranking <- ind_ranking[order(-ind_ranking$ind_prop), ]
  ind_ranking$ind_rank <- 1:nrow(ind_ranking)
  prs_ranking_summary <- merge(pop_ranking[, c("PRS", "prop", "pop_rank")],
                                ind_ranking[, c("PRS", "ind_prop", "ind_rank")], by = "PRS")
  ranking <- prs_ranking_summary[, c("PRS", "pop_rank", "ind_rank")]
  colnames(ranking) <- c("PRS", "Population", "Individual")
}

ranking <- ranking %>%
  mutate(slope_dir = case_when(
    Individual < Population ~ "Positive",
    Individual > Population ~ "Negative",
    TRUE ~ "Constant"
  ))

ranking_long <- ranking %>%
  pivot_longer(cols = c("Population", "Individual"),
               names_to = "Rank_Type", values_to = "Rank")

p_slope <- ggplot(ranking_long, aes(x = Rank_Type, y = Rank, group = PRS)) +
  geom_line(aes(color = slope_dir), linewidth = 1.2) +
  geom_point(aes(color = slope_dir), size = 4) +
  scale_y_reverse() +
  scale_color_manual(values = c(
    "Positive" = "red",
    "Negative" = "blue",
    "Constant" = "gray"
  )) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  geom_text(
    data = ranking_long %>% filter(Rank_Type == "Individual"),
    aes(label = PRS, x = 1),
    hjust = 1.1, size = 3
  ) +
  geom_text(
    data = ranking_long %>% filter(Rank_Type == "Population"),
    aes(label = PRS, x = 2),
    hjust = -0.1, size = 3
  ) +
  coord_cartesian(clip = 'off')

ggsave(file.path(OUTPUT_DIR, "7_plots", "myplot_slope.png"),
       plot = p_slope, width = 8, height = 12, units = "in", dpi = 300, bg = "white")

report_step("Step 7b")
} else {
  cat("\n========== STEP 7: Skipped ==========\n")
}

###############################################################################
total_elapsed <- (proc.time() - pipeline_start)["elapsed"]
peak_mem_mb <- as.numeric(gc(verbose = FALSE)[2, 2])
n_output_files <- length(list.files(OUTPUT_DIR, recursive = TRUE, pattern = "\\.(csv|png|txt|xlsx)$"))

cat("\n========== PIPELINE COMPLETE ==========\n")
cat("All outputs saved to:", OUTPUT_DIR, "\n\n")
cat("--- Run Summary ---\n")
cat(sprintf("  Total runtime:   %.1f sec (%.1f min)\n", total_elapsed, total_elapsed / 60))
cat(sprintf("  Peak memory:     %.0f MB\n", peak_mem_mb))
cat(sprintf("  Output files:    %d\n", n_output_files))
cat(sprintf("  R version:       %s\n", R.version.string))
cat(sprintf("  Platform:        %s\n", R.version$platform))
cat(sprintf("  Random seed:     42\n"))
cat(sprintf("  CPUs available:  %d\n", parallel::detectCores()))
cat("-------------------\n")

# Write summary to file
summary_df <- data.frame(
  metric = c("total_runtime_sec", "peak_memory_mb", "output_files",
             "r_version", "platform", "random_seed", "cpus"),
  value = c(round(total_elapsed, 1), round(peak_mem_mb), n_output_files,
            R.version.string, R.version$platform, "42",
            parallel::detectCores())
)
write.csv(summary_df, file.path(OUTPUT_DIR, "run_summary.csv"), row.names = FALSE)
