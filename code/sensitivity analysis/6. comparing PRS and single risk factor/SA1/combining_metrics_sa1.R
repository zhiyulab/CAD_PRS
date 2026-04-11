######################################################################
library(data.table)
library(readr)

pheno<-fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt") 
pheno = data.frame(pheno)
pheno$IDs_toRemove_SampleQCv3 = ifelse((!(pheno$Submitted_Gender == pheno$Inferred_Gender) | pheno$Non_Consented== 1),1,0) 
pheno = pheno[which(pheno$IDs_toRemove_SampleQCv3==0),]

#get lpa and crp colulmns
lpa<-pheno[,c(1,188)]
colnames(lpa)<-c("id","lpa")
lpa <- lpa[complete.cases(lpa), ]


pce_zy_2 <- read_delim("/medpop/esp2/sliang/pce_zy_2.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
PCE <- pce_zy_2[, c(1,12)]
colnames(PCE)<-c("id","pce")
PCE <- PCE[complete.cases(PCE), ]

#standardize
UKB_PRS <- read_csv("/home/unix/sliang/sliang/age/UKB_PRS_58.csv")
UKB_PRS <- UKB_PRS[complete.cases(UKB_PRS), ]
UKB_PRS_lpa <- merge(UKB_PRS, lpa, by = "id")
UKB_PRS_lpa_pce <- merge(UKB_PRS_lpa, PCE, by = "id")


UKB_PRS_lpa_pce$lpa_standardized <- as.numeric(scale(UKB_PRS_lpa_pce$lpa))



for (i in 1:58) {
  col_name <- paste0("df", i)
  new_col  <- paste0(col_name, "_lpa")
  
  UKB_PRS_lpa_pce[[new_col]] <- UKB_PRS_lpa_pce[[col_name]] + UKB_PRS_lpa_pce$lpa_standardized
}

for (i in 1:58) {
  col_name <- paste0("df", i)
  new_col  <- paste0(col_name, "_pce")
  
  UKB_PRS_lpa_pce[[new_col]] <- UKB_PRS_lpa_pce[[col_name]] + UKB_PRS_lpa_pce$pce
}


for (i in 1:58) {
  col_name <- paste0("df", i, "_lpa")
  
  cutoff <- quantile(UKB_PRS_lpa_pce[[col_name]], 0.9, na.rm = TRUE)
  new_col_name <- paste0(col_name, "_cutoff")
  
  UKB_PRS_lpa_pce[[new_col_name]] <- ifelse(
    UKB_PRS_lpa_pce[[col_name]] >= cutoff, 1, 0
  )
}

for (i in 1:58) {
  col_name <- paste0("df", i, "_pce")
  
  cutoff <- quantile(UKB_PRS_lpa_pce[[col_name]], 0.9, na.rm = TRUE)
  new_col_name <- paste0(col_name, "_cutoff")
  
  UKB_PRS_lpa_pce[[new_col_name]] <- ifelse(
    UKB_PRS_lpa_pce[[col_name]] >= cutoff, 1, 0
  )
}


cohort_p_55 <- read.csv("/home/unix/sliang/sliang/age/cohort_p_sex_specific.csv")
UKB_PRS_lpa_pce_p<-UKB_PRS_lpa_pce[UKB_PRS_lpa_pce$id %in% cohort_p_55$id,]
fwrite(UKB_PRS_lpa_pce_p, "UKB_PRS_lpa_pce_p.csv")


###############################################################

library(boot)

# Select the column
x <- UKB_PRS_lpa_pce_p$df53_lpa_cutoff
x <- x[!is.na(x)]  # remove NAs

# Run bootstrap
boot_out <- boot(
  data = x,
  statistic = function(x, indices) {
    mean(x[indices] == 1)  # proportion of 1s
  },
  R = 1000  # number of bootstrap resamples
)

# Compute percentile confidence interval
ci <- boot.ci(boot_out, type = "perc")

# Collect results in a data frame
results_prs_lpa <- data.frame(
  variable = "prs_lpa_cutoff",
  prop = mean(x == 1),
  lower_95 = ci$percent[4],
  upper_95 = ci$percent[5]
)

fwrite(results_prs_lpa,"results_prs_lpa.csv")
###############################################################

library(boot)

# Select the column
x <- UKB_PRS_lpa_pce_p$df53_pce_cutoff
x <- x[!is.na(x)]  # remove NAs

# Run bootstrap
boot_out <- boot(
  data = x,
  statistic = function(x, indices) {
    mean(x[indices] == 1)  # proportion of 1s
  },
  R = 1000  # number of bootstrap resamples
)

# Compute percentile confidence interval
ci <- boot.ci(boot_out, type = "perc")

# Collect results in a data frame
results_prs_pce <- data.frame(
  variable = "prs_pce_cutoff",
  prop = mean(x == 1),
  lower_95 = ci$percent[4],
  upper_95 = ci$percent[5]
)

fwrite(results_prs_pce,"results_prs_pce.csv")