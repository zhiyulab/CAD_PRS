#read 58 PRS files 
# Step 1: Directory
data_dir <- "/medpop/esp2/yang/project/cli_trail/CAD/UKB/resultsum_SQ"

# Step 2: List all matching files
all_files <- list.files(path = data_dir, pattern = "^PC4\\.(.*)\\.adjNormPRS\\.txt\\.gz$", full.names = TRUE)

# Step 3: Exclude specific PGS IDs
excluded_ids <- c("PGS000116", "PGS004899", "PGS004888", "PGS004237")

# Step 4: Filter files by excluding those with unwanted IDs
filtered_files <- all_files[!grepl(paste(excluded_ids, collapse = "|"), all_files)]

# Step 5: sort to maintain consistent ordering
filtered_files <- sort(filtered_files)

# Step 6: Read files into a list
data_list <- lapply(filtered_files, function(f) read.table(f, header = TRUE, sep = "\t"))

# Step 7: Assign to df1, df2, ...
for (i in seq_along(data_list)) {
  assign(paste0("df", i), data_list[[i]])
}
#################################################################################
#extract only 4th column from df1 and df58 and merge into one file
# Initialize an empty list to hold the 4th columns
prs_list <- list()

# Loop through df1 to df58
for (i in 1:58) {
  df <- get(paste0("df", i))
  prs_list[[i]] <- df[[4]]  # 4th column
}

# Combine all into one data frame (columns side-by-side)
UKB_PRS_58 <- as.data.frame(prs_list)

# Rename columns as df1, df2, ...
colnames(UKB_PRS_58) <- paste0("df", 1:58)

UKB_PRS_58$id<-df1$IID
UKB_PRS_58 <- UKB_PRS_58[, c(ncol(UKB_PRS_58), 1:(ncol(UKB_PRS_58) - 1))]

#################################################################################
#for each score, individuals with top 10% score out of whole UKBB is flagged cutoff=1
for (i in 1:58) {
  col_name <- paste0("df", i)
  cutoff <- quantile(UKB_PRS_58[[col_name]], 0.9, na.rm = TRUE)
  new_col_name <- paste0(col_name, "_cutoff")
  
  UKB_PRS_58[[new_col_name]] <- ifelse(UKB_PRS_58[[col_name]] >= cutoff, 1, 0)
}
#################################################################################
#only include scores and cutoff status for SMuRFLess early onset individuals
library(readr)
library(data.table)
cohort_p <- read.csv("/home/unix/sliang/not_cohort_p.csv")
cohort_p<-cohort_p[,c(1,10)]
cohort_p<-merge(cohort_p,UKB_PRS_58,by="id",all.x=TRUE)
cohort_p_other<-cohort_p[(cohort_p$in_white_British_ancestry_subset==0),]
cohort_p_white<-cohort_p[(cohort_p$in_white_British_ancestry_subset==1),]

##############################################testing bootstrap
library(boot)

n_boot <- 1000
results <- lapply(1:58, function(i) {
  col_name <- paste0("df", i, "_cutoff")
  x <- cohort_p[[col_name]]
  x <- x[!is.na(x)]
  
  boot_out <- boot(
    data = x,
    statistic = function(x, indices) {
      mean(x[indices] == 1)
    },
    R = 1000
  )
  
  ci <- boot.ci(boot_out, type = "perc")
  
  data.frame(
    PRS = paste0("df", i),
    prop = mean(x == 1),
    lower_95 = ci$percent[4],
    upper_95 = ci$percent[5]
  )
})
boot_summary_df_55 <- do.call(rbind, results)
colnames(boot_summary_df_55) <- c("PRS", "prop", "low_ci", "high_ci")
library(data.table)
fwrite(boot_summary_df_55, "boot_summary_df_55.csv")

##############################################testing bootstrap
library(boot)

n_boot <- 1000
results <- lapply(1:58, function(i) {
  col_name <- paste0("df", i, "_cutoff")
  x <- cohort_p_white[[col_name]]
  x <- x[!is.na(x)]
  
  boot_out <- boot(
    data = x,
    statistic = function(x, indices) {
      mean(x[indices] == 1)
    },
    R = 1000
  )
  
  ci <- boot.ci(boot_out, type = "perc")
  
  data.frame(
    PRS = paste0("df", i),
    prop = mean(x == 1),
    lower_95 = ci$percent[4],
    upper_95 = ci$percent[5]
  )
})
boot_summary_df_55 <- do.call(rbind, results)
colnames(boot_summary_df_55) <- c("PRS", "prop", "low_ci", "high_ci")
library(data.table)
fwrite(boot_summary_df_55, "boot_summary_df_55_white.csv")

##############################################testing bootstrap
library(boot)

n_boot <- 1000
results <- lapply(1:58, function(i) {
  col_name <- paste0("df", i, "_cutoff")
  x <- cohort_p_other[[col_name]]
  x <- x[!is.na(x)]
  
  boot_out <- boot(
    data = x,
    statistic = function(x, indices) {
      mean(x[indices] == 1)
    },
    R = 1000
  )
  
  ci <- boot.ci(boot_out, type = "perc")
  
  data.frame(
    PRS = paste0("df", i),
    prop = mean(x == 1),
    lower_95 = ci$percent[4],
    upper_95 = ci$percent[5]
  )
})
boot_summary_df_55 <- do.call(rbind, results)
colnames(boot_summary_df_55) <- c("PRS", "prop", "low_ci", "high_ci")
library(data.table)
fwrite(boot_summary_df_55, "boot_summary_df_55_other.csv")
stop()
######################################################################
