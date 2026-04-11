# Load libraries
if (!require("data.table")) install.packages("data.table", repos="https://cloud.r-project.org")
library(data.table)

# === Set working directory and input paths ===
output_dir <- "/medpop/esp2/sliang"

# === Load phenotype data ===
pheno<-fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt") 

# === Quality control: remove mismatched gender or non-consented samples ===
pheno$IDs_toRemove_SampleQCv3 <- ifelse(
  (!(pheno$Submitted_Gender == pheno$Inferred_Gender) | pheno$Non_Consented == 1),
  1, 0
)
pheno <- pheno[pheno$IDs_toRemove_SampleQCv3 == 0, ]

# === Load second phenotype dataset for intermediate CAD ===
X20250615_pheno_noGP <- fread("/medpop/esp2/lli/UKB_pheno_curation/Results/20250615_pheno_noGP.tsv.gz")
library(readr)
pce_zy_2 <- read_delim("/medpop/esp2/sliang/pce_zy_2.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
PCE <- pce_zy_2[, c(1,12)]

# === Select individuals with intermediate CAD ===
ID_intermediate <- subset(X20250615_pheno_noGP, has_disease_Coronary_Artery_Disease_INTERMEDIATE == 1)
nrow(ID_intermediate)
subset_data_intermediate <- pheno[pheno$id %in% ID_intermediate$sample_id, ]
nrow(subset_data_intermediate)
# === Select key phenotype columns (Smoking, SBP, DBP, HbA1c, Choletserol, LDL, Medication) and remove missing values ===
data_intermediate <- subset_data_intermediate[, c(1,12,14,15,23,22,74,184,178,187,2226,2227,2228)]


data_intermediate<- merge(data_intermediate,PCE, by = "id", all.x = TRUE)
nrow(data_intermediate)

data_intermediate <- data_intermediate[complete.cases(data_intermediate), ]
nrow(data_intermediate)

colnames(data_intermediate)[ncol(data_intermediate)] <- "PCE"

# === Convert HbA1c to NGSP units ===
data_intermediate$HbA1c <- data_intermediate$Glycated.haemoglobin..HbA1c. * 0.09148 + 2.152

# === Remove individuals with SMuRFs (conventional risk factors defined by lab test results and medication status) ===
data_intermediate$IDs_toRemove <- ifelse(
  data_intermediate$HbA1c >= 6.5 |
    data_intermediate$PCE >= 0.05 |
    data_intermediate$LDL.direct >= 4.9 |
    data_intermediate$Diabetes_Med == 1 |
    data_intermediate$LipidLowering_Med == 1,
  1, 0
)
cohort_p <- subset(data_intermediate, IDs_toRemove == 0)


# === Calculate onset age and restrict to early onset <55 ===
time <- X20250615_pheno_noGP[X20250615_pheno_noGP$sample_id %in% cohort_p$id, ]
time <- time[, c(1,178,179,180,181)]
colnames(time)[1] <- "id"

cohort_p_all <- merge(time, cohort_p, by = "id", all.x = TRUE)
cohort_p_all$onset_age <- cohort_p_all$age + cohort_p_all$followup_Coronary_Artery_Disease_INTERMEDIATE
cohort_p_sex_specific <- subset(
  cohort_p_all,
  (Sex == "Male"   & onset_age < 55) |
    (Sex == "Female" & onset_age < 60)
)


fwrite(cohort_p_sex_specific, "cohort_p_sex_specific.csv")
