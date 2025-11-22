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
output_dir <- "/medpop/esp2/sliang"
fwrite(UKB_PRS_58, file.path(output_dir, "UKB_PRS.csv"))
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
SMuRFless_55 <- read.csv("/medpop/esp2/sliang/SMuRFless_CAD_under55.csv")
race<-SMuRFless_55[,c("id","in_white_British_ancestry_subset")]
SMuRFless_55<-merge(SMuRFless_55,UKB_PRS_58,by="id",all.x=TRUE)
SMuRFless_55<-SMuRFless_55[,c(1,2330:2445)]
SMuRFless_55<-merge(SMuRFless_55,race, by="id",all.x=TRUE)
SMuRFless_55_other<-SMuRFless_55[(SMuRFless_55$in_white_British_ancestry_subset==0),]
SMuRFless_55_white<-SMuRFless_55[(SMuRFless_55$in_white_British_ancestry_subset==1),]

#################################################################################
#calculate percentage of SMuRFless early onset flagged as cutoff=1
#(number of cutoff==1) / (size of SMuRFless early onset group - number of no score)

# Initialize vectors to store results
na_counts <- numeric(58)
one_counts <- numeric(58)

# Loop through columns df1_5 to df58_5
for (i in 1:58) {
  col_name <- paste0("df", i, "_cutoff")
  column_data <- SMuRFless_55[[col_name]]
  
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

summary_df_55$PRS <- gsub("^df", "", summary_df_55$PRS)
summary_df_55$PRS<-as.numeric(summary_df_55$PRS)
summary_df_55 <- summary_df_55[order(summary_df_55$PRS), ]
rownames(summary_df_55) <- as.character(1:nrow(summary_df_55))
summary_df_55$all_percentage<-summary_df_55$SMuRFless_cutoff1/(nrow(SMuRFless_55)-summary_df_55$SMuRFless_NA)
#################################################################################
#repeat for white individuals out of SMuRFless early onset

# Initialize vectors to store results
na_counts <- numeric(58)
one_counts <- numeric(58)

# Loop through columns df1_10 to df58_10
for (i in 1:58) {
  col_name <- paste0("df", i, "_cutoff")
  column_data <- SMuRFless_55_white[[col_name]]
  
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
summary_df_55_white$white_percentage<-summary_df_55_white$SMuRFless_55_white_cutoff1/(nrow(SMuRFless_55_white)-summary_df_55$SMuRFless_55_white_NA)

summary_df_55<-merge(summary_df_55, summary_df_55_white,by = "PRS", all.x = TRUE)

#################################################################################
#repeat for other ancestry individuals out of SMuRFless early onset
# Initialize vectors to store results
na_counts <- numeric(58)
one_counts <- numeric(58)

# Loop through columns df1_10 to df58_10
for (i in 1:58) {
  col_name <- paste0("df", i, "_cutoff")
  column_data <- SMuRFless_55_other[[col_name]]
  
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
summary_df_55_other$other_percentage<-summary_df_55_other$SMuRFless_55_other_cutoff1/((nrow(SMuRFless_55_other)-summary_df_55$SMuRFless_55_other_NA))

summary_df_55<-merge(summary_df_55, summary_df_55_other,by = "PRS", all.x = TRUE)
library(data.table)
output_dir <- "/medpop/esp2/sliang"
fwrite(summary_df_55, file.path(output_dir, "summary_df_55_final.csv"))


#calculate number of pick ups (number of cutoff=1) for each individual
#only include cutoff status for each individual
SMuRFless_55_sum<-SMuRFless_55[,c(1,60:117)]
#sum is the number of 1s each individual has
SMuRFless_55_sum$sum <- rowSums(SMuRFless_55_sum[, !names(SMuRFless_55_sum) %in% "id"], na.rm = TRUE)
fwrite(SMuRFless_55_sum, file.path(output_dir, "SMuRFless_55_sum.csv"))

