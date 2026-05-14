# Step 1: Define the directory where the files are stored
data_dir <- "~/Partners HealthCare Dropbox/Shengxin Liang/sliang/PRS_ys/resultsum_SQ"

# Step 2: List all matching files
all_files <- list.files(path = data_dir, pattern = "^PC4\\.(.*)\\.adjNormPRS\\.txt\\.gz$", full.names = TRUE)

# Step 3: Exclude specific PGS IDs
excluded_ids <- c("PGS000116", "PGS004899", "PGS004888", "PGS004237")

# Step 4: Filter files by excluding those with unwanted IDs
filtered_files <- all_files[!grepl(paste(excluded_ids, collapse = "|"), all_files)]

# Optional: sort to maintain consistent ordering
filtered_files <- sort(filtered_files)

# Step 5: Read files into a list
data_list <- lapply(filtered_files, function(f) read.table(f, header = TRUE, sep = "\t"))

# Step 6: Assign to df1, df2, ...
for (i in seq_along(data_list)) {
  assign(paste0("df", i), data_list[[i]])
}
#################################################################################
# Initialize an empty list to hold the 4th columns
prs_list <- list()

# Loop through df1 to df58
for (i in 1:58) {
  df <- get(paste0("df", i))
  prs_list[[i]] <- df[[4]]  # 4th column
}

# Combine all into one data frame (columns side-by-side)
MGB_PRS_58 <- as.data.frame(prs_list)

# Optional: Rename columns as df1, df2, ...
colnames(MGB_PRS_58) <- paste0("df", 1:58)

MGB_PRS_58$IID<-df$IID
MGB_PRS_58 <- MGB_PRS_58[, c(ncol(MGB_PRS_58), 1:(ncol(MGB_PRS_58) - 1))]
MGB_PRS_58$IID <- sub(".*-", "", MGB_PRS_58$IID)
colnames(MGBB_Genomic_Data_Linker_2022_06_09)[1] <- "IID"
MGB_PRS_58<-merge(MGB_PRS_58, MGBB_Genomic_Data_Linker_2022_06_09, by="IID",all.x=T)
MGB_PRS_58$IID<-NULL
MGB_PRS_58$EPIC_PMRN<-NULL
MGB_PRS_58$MGH_MRN<-NULL
MGB_PRS_58 <- MGB_PRS_58[, c(ncol(MGB_PRS_58), 1:(ncol(MGB_PRS_58) - 1))]

rm(list = paste0("df", 1:58))
#################################################################################
# Let's assume UKB_PRS_58 has columns named df1 to df58
for (i in 1:58) {
  col_name <- paste0("df", i)
  cutoff <- quantile(MGB_PRS_58[[col_name]], 0.9, na.rm = TRUE)
  new_col_name <- paste0(col_name, "_cutoff")
  
  MGB_PRS_58[[new_col_name]] <- ifelse(MGB_PRS_58[[col_name]] >= cutoff, 1, 0)
}
#################################################################################

SMuRFless_55<-smurfless_early
MGB_55_PRS<-merge(SMuRFless_55,MGB_PRS_58,by="EMPI",all.x=T)


MGB_55_female <- MGB_55_PRS[MGB_55_PRS$Sex == "Female",]
MGB_55_Male <- MGB_55_PRS[MGB_55_PRS$Sex == "Male",]

#################################################################################
#11.3% of SMuRFless group are in top 5%?
#(no.cutoff==1 in SMuRFless) / (no.SMuRFless - no.NA in SMuRFless)
#calculate number of n/a and cutoff=1 for SMuRFless_55 all ancestry
# Initialize vectors to store results
na_counts <- numeric(58)
one_counts <- numeric(58)

# Loop through columns df1_5 to df58_5
for (i in 1:58) {
  col_name <- paste0("df", i, "_cutoff")
  column_data <- MGB_55_PRS[[col_name]]
  
  na_counts[i] <- sum(is.na(column_data))
  one_counts[i] <- sum(column_data == 1, na.rm = TRUE)
}

# Create summary data frame
summary_df_55 <- data.frame(
  PRS = paste0("df", 1:58),
  NA_count = na_counts,
  One_count = one_counts
)
colnames(summary_df_55)<-c("PRS","SMuRFless_NA","SMuRFless_cutoff1")

##############################################testing bootstrap
library(boot)
 
n_boot <- 1000
results <- lapply(1:58, function(i) {
  col_name <- paste0("df", i, "_cutoff")
  x <- MGB_55_Male[[col_name]]
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
fwrite(boot_summary_df_55, "boot_summary_df_55_male.csv")

#########################################################

summary_df_55$PRS <- gsub("^df", "", summary_df_55$PRS)
summary_df_55$PRS<-as.numeric(summary_df_55$PRS)
summary_df_55 <- summary_df_55[order(summary_df_55$PRS), ]
rownames(summary_df_55) <- as.character(1:nrow(summary_df_55))
summary_df_55$all_percentage<-summary_df_55$SMuRFless_cutoff1/(661-160)
#################################################################################
#repeat for SMuRFless_55_white

# Initialize vectors to store results
na_counts <- numeric(58)
one_counts <- numeric(58)

# Loop through columns df1_10 to df58_10
for (i in 1:58) {
  col_name <- paste0("df", i, "_cutoff")
  column_data <- MGB_55_white[[col_name]]
  
  na_counts[i] <- sum(is.na(column_data))
  one_counts[i] <- sum(column_data == 1, na.rm = TRUE)
}

# Create summary data frame
summary_df_55_white <- data.frame(
  PRS = paste0("df", 1:58),
  NA_count = na_counts,
  One_count = one_counts
)
colnames(summary_df_55_white)<-c("PRS","SMuRFless_55_white_NA","SMuRFless_55_white_cutoff1")

summary_df_55_white$PRS <- gsub("^df", "", summary_df_55_white$PRS)
summary_df_55_white$PRS<-as.numeric(summary_df_55_white$PRS)
summary_df_55_white <- summary_df_55_white[order(summary_df_55_white$PRS), ]
rownames(summary_df_55_white) <- as.character(1:nrow(summary_df_55_white))
summary_df_55_white$white_percentage<-summary_df_55_white$SMuRFless_55_white_cutoff1/(558-142)

summary_df_55<-merge(summary_df_55, summary_df_55_white,by = "PRS", all.x = TRUE)

#################################################################################
#repeat for SMuRFless_55_other
# Initialize vectors to store results
na_counts <- numeric(58)
one_counts <- numeric(58)

# Loop through columns df1_10 to df58_10
for (i in 1:58) {
  col_name <- paste0("df", i, "_cutoff")
  column_data <- MGB_55_other[[col_name]]
  
  na_counts[i] <- sum(is.na(column_data))
  one_counts[i] <- sum(column_data == 1, na.rm = TRUE)
}

# Create summary data frame
summary_df_55_other <- data.frame(
  PRS = paste0("df", 1:58),
  NA_count = na_counts,
  One_count = one_counts
)
colnames(summary_df_55_other)<-c("PRS","SMuRFless_55_other_NA","SMuRFless_55_other_cutoff1")

summary_df_55_other$PRS <- gsub("^df", "", summary_df_55_other$PRS)
summary_df_55_other$PRS<-as.numeric(summary_df_55_other$PRS)
summary_df_55_other <- summary_df_55_other[order(summary_df_55_other$PRS), ]
rownames(summary_df_55_other) <- as.character(1:nrow(summary_df_55_other))
summary_df_55_other$other_percentage<-summary_df_55_other$SMuRFless_55_other_cutoff1/(73-18)

summary_df_55<-merge(summary_df_55, summary_df_55_other,by = "PRS", all.x = TRUE)
write.csv(summary_df_55, "summary_df_55_final.csv", row.names = FALSE)

