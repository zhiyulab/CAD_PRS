
# Step 1: Directory
data_dir <- "/medpop/esp2/yang/project/for_others/sliang/3scores/resultsum_SQ"

# Step 2: List all matching files
all_files <- list.files(
  path = data_dir,
  pattern = "^PC4\\.(.*)\\.adjNormPRS\\.txt$",
  full.names = TRUE
)

# Step 3: Sort to maintain consistent ordering
all_files <- sort(all_files)

# Step 4: Read files into a list
data_list <- lapply(all_files, function(f) read.table(f, header = TRUE, sep = "\t"))

# Step 5: Assign to df1, df2, ...
for (i in seq_along(data_list)) {
  assign(paste0("df", i), data_list[[i]])
}
#################################################################################
#extract only 4th column from df1 and df58 and merge into one file
# Initialize an empty list to hold the 4th columns
prs_list <- list()

# Loop through df1 to df58
for (i in 1:3) {
  df <- get(paste0("df", i))
  prs_list[[i]] <- df[[4]]  # 4th column
}

# Combine all into one data frame (columns side-by-side)
lipid_removed_prs <- as.data.frame(prs_list)

# Rename columns as df1, df2, ...
colnames(lipid_removed_prs) <- paste0("df", 1:3)

lipid_removed_prs$id<-df1$IID
lipid_removed_prs <- lipid_removed_prs[, c(ncol(lipid_removed_prs), 1:(ncol(lipid_removed_prs) - 1))]

#################################################################################
#for each score, individuals with top 10% score out of whole UKBB is flagged cutoff=1
for (i in 1:3) {
  col_name <- paste0("df", i)
  cutoff <- quantile(lipid_removed_prs[[col_name]], 0.9, na.rm = TRUE)
  new_col_name <- paste0(col_name, "_cutoff")
  
  lipid_removed_prs[[new_col_name]] <- ifelse(lipid_removed_prs[[col_name]] >= cutoff, 1, 0)
}
library(data.table)
fwrite(lipid_removed_prs, "lipid_removed_prs.csv")

#################################################################################
#only include scores and cutoff status for SMuRFLess early onset individuals
library(readr)
library(data.table)
cohort_p <- read.csv("/home/unix/sliang/cohort_p_55.csv")
cohort_p<-cohort_p[,c(1,10)]
cohort_p<-merge(cohort_p,lipid_removed_prs,by="id",all.x=TRUE)
cohort_p_other<-cohort_p[(cohort_p$in_white_British_ancestry_subset==0),]
cohort_p_white<-cohort_p[(cohort_p$in_white_British_ancestry_subset==1),]
fwrite(cohort_p, "cohort_p.csv")
fwrite(cohort_p_other, "cohort_p_other.csv")
fwrite(cohort_p_white, "cohort_p_white.csv")

##############################################testing bootstrap
library(boot)

n_boot <- 1000
results <- lapply(1:3, function(i) {
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
results <- lapply(1:3, function(i) {
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
results <- lapply(1:3, function(i) {
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
