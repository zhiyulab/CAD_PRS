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

# === Select individuals with intermediate CAD ===
ID_intermediate <- subset(X20250615_pheno_noGP, has_disease_Coronary_Artery_Disease_INTERMEDIATE == 1)
subset_data_intermediate <- pheno[pheno$id %in% ID_intermediate$sample_id, ]

# === Select key phenotype columns (Smoking, SBP, DBP, HbA1c, Choletserol, LDL, Medication) and remove missing values ===
data_intermediate <- subset_data_intermediate[, c(1,14,15,23,22,184,178,187,2226,2227,2228)]
data_intermediate <- data_intermediate[complete.cases(data_intermediate), ]

# === Convert HbA1c to NGSP units ===
data_intermediate$HbA1c <- data_intermediate$Glycated.haemoglobin..HbA1c. * 0.09148 + 2.152

# === Remove individuals with SMuRFs (conventional risk factors defined by lab test results and medication status) ===
data_intermediate$IDs_toRemove <- ifelse(
  data_intermediate$SmokingStatus == 2 |
    data_intermediate$SBP >= 140 |
    data_intermediate$DBP >= 90 |
    data_intermediate$HbA1c >= 6.5 |
    data_intermediate$Cholesterol > 5.5 |
    data_intermediate$LDL.direct > 3.5 |
    data_intermediate$Diabetes_Med == 1 |
    data_intermediate$Antihypertensive_Med == 1 |
    data_intermediate$LipidLowering_Med == 1,
  1, 0
)
SMuRFless_intermediate <- subset(data_intermediate, IDs_toRemove == 0)

# === Identify individuals whose CAD occurred before other chronic conditions using follw up time of diagnosis===
time <- X20250615_pheno_noGP[X20250615_pheno_noGP$sample_id %in% SMuRFless_intermediate$id, ]
time <- time[, c(1,178,179,180,181,214,215,216,217,270,271,272,273,266,267,268,269)]
colnames(time)[1] <- "id"

# Step 1: remove if CAD occurred after other conditions
time$IDs_toRemove <- ifelse(
  time$followup_Coronary_Artery_Disease_INTERMEDIATE >= time$followup_Diabetes_All |
    time$followup_Coronary_Artery_Disease_INTERMEDIATE >= time$followup_Hypercholesterolemia |
    time$followup_Coronary_Artery_Disease_INTERMEDIATE >= time$followup_Hypertension,
  1, 0
)
time_kept <- subset(time, IDs_toRemove == 0)

# Step 2: remove if time difference < 0.5 years
time_kept$IDs_toRemove <- ifelse(
  (time_kept$followup_Diabetes_All - time_kept$followup_Coronary_Artery_Disease_INTERMEDIATE < 0.5) |
    (time_kept$followup_Hypercholesterolemia - time_kept$followup_Coronary_Artery_Disease_INTERMEDIATE < 0.5) |
    (time_kept$followup_Hypertension - time_kept$followup_Coronary_Artery_Disease_INTERMEDIATE < 0.5),
  1, 0
)
time_final <- subset(time_kept, IDs_toRemove == 0)

# === Calculate onset age and restrict to early onset <55 ===
SMuRFless_all <- merge(time_final, subset_data_intermediate, by = "id", all.x = TRUE)
SMuRFless_all$onset_age <- SMuRFless_all$age + SMuRFless_all$followup_Coronary_Artery_Disease_INTERMEDIATE
SMuRFless_55 <- subset(SMuRFless_all, onset_age < 55)

# === Save outputs ===
fwrite(SMuRFless_55, file.path(output_dir, "SMuRFless_CAD_under55.csv"))
fwrite(SMuRFless_all,file.path(output_dir, "SMuRFless_CAD.csv"))
complete_data_intermediate<-subset_data_intermediate[subset_data_intermediate$id %in% data_intermediate$id, ]
SMuRF <- complete_data_intermediate[!complete_data_intermediate$id %in% SMuRFless_all$id, ]
fwrite(SMuRF,file.path(output_dir, "SMuRF.csv"))