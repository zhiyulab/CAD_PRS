if (!require("tableone")) install.packages("tableone", repos="https://cloud.r-project.org")
library(tableone)
if (!require("readr")) install.packages("readr", repos = "https://cloud.r-project.org")
library(readr)
library(data.table)
LE8score_composite_ <- read_csv("/medpop/esp2/yxliu/LE8score_calculation/updated_LE8score_wComposite.txt")
pce_zy_2 <- read_delim("/medpop/esp2/sliang/pce_zy_2.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
PCE <- pce_zy_2[, c(1,12)]
family<-read_csv("/home/unix/sliang/sliang/family/family_history_final.csv")
fh<-read_delim("/medpop/esp2/sliang/ukb_exome_450k_fh.carrier.txt", 
               delim = "\t", escape_double = FALSE, 
               trim_ws = TRUE)
fh<-fh[,c(1,3)]
colnames(fh)[1] <- "id"

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

########################################################
#get columns with variables of interest to present in table one
sa1<-read_csv("/home/unix/sliang/sliang/age/cohort_p_sex_specific.csv")
sa1<-sa1[,"id"]
tableone_sa1<-merge(sa1, pheno, by = "id", all.x = TRUE)
tableone_sa1<-tableone_sa1[, c(1,12,14,74,180,2309)]
tableone_sa1 <- merge(tableone_sa1, PCE, by = "id", all.x = TRUE)
tableone_sa1 <- merge(tableone_sa1, LE8score_composite_, by = "id", all.x = TRUE)
tableone_sa1 <- merge(tableone_sa1, family, by = "id", all.x = TRUE)
tableone_sa1 <- merge(tableone_sa1, fh, by = "id", all.x = TRUE)

#remove repeated columns
tableone_sa1$Sex.y<-NULL
tableone_sa1$age.y<-NULL

colnames(tableone_sa1) 

colnames(tableone_sa1) <- c("id", "Sex", "Age", "White British Ancestry", "C Reactive Protein","Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                            "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")

tab4 <- CreateTableOne(vars = myVars, data = tableone_sa1, factorVars = catVars)
table1_sa1<-print(tab4,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_sa1, file = "table1_sa1.csv")
########################################################
#get columns with variables of interest to present in table one
sa2<-read_csv("/home/unix/sliang/sliang/family/family_history_both_CAD.csv")
sa2<-sa2[,"id"]
tableone_sa2<-merge(sa2, pheno, by = "id", all.x = TRUE)
tableone_sa2<-tableone_sa2[, c(1,12,14,74,180,2309)]
tableone_sa2 <- merge(tableone_sa2, PCE, by = "id", all.x = TRUE)
tableone_sa2 <- merge(tableone_sa2, LE8score_composite_, by = "id", all.x = TRUE)
tableone_sa2 <- merge(tableone_sa2, family, by = "id", all.x = TRUE)
tableone_sa2 <- merge(tableone_sa2, fh, by = "id", all.x = TRUE)

#remove repeated columns
tableone_sa2$Sex.y<-NULL
tableone_sa2$age.y<-NULL

colnames(tableone_sa2) <- c("id", "Sex", "Age", "White British Ancestry", "C Reactive Protein","Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                            "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")

tab5 <- CreateTableOne(vars = myVars, data = tableone_sa2, factorVars = catVars)
table1_sa2<-print(tab5,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_sa2, file = "table1_sa2.csv")
########################################################
#get columns with variables of interest to present in table one
family_cohort <- read_csv("/home/unix/sliang/sliang/family/family_history_both_CAD.csv")
cohort_p<-read_csv("/home/unix/sliang/sliang/family/cohort_p_55.csv")
sa3<-cohort_p[cohort_p$id %in% family_cohort$id,]
sa3<-sa3[,"id"]
tableone_sa3<-merge(sa3, pheno, by = "id", all.x = TRUE)
tableone_sa3<-tableone_sa3[, c(1,12,14,74,180,2309)]
tableone_sa3 <- merge(tableone_sa3, PCE, by = "id", all.x = TRUE)
tableone_sa3 <- merge(tableone_sa3, LE8score_composite_, by = "id", all.x = TRUE)
tableone_sa3 <- merge(tableone_sa3, family, by = "id", all.x = TRUE)
tableone_sa3 <- merge(tableone_sa3, fh, by = "id", all.x = TRUE)

#remove repeated columns
tableone_sa3$Sex.y<-NULL
tableone_sa3$age.y<-NULL

colnames(tableone_sa3) <- c("id", "Sex", "Age", "White British Ancestry", "C Reactive Protein","Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                            "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")

tab6 <- CreateTableOne(vars = myVars, data = tableone_sa3, factorVars = catVars)
table1_sa3<-print(tab6,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_sa3, file = "table1_sa3.csv")
########################################################
#get columns with variables of interest to present in table one
sa4<-read_csv("/home/unix/sliang/sliang/smurfless/SMuRFless_CAD_under55.csv")
sa4<-sa4[,"id"]
tableone_sa4<-merge(sa4, pheno, by = "id", all.x = TRUE)
tableone_sa4<-tableone_sa4[, c(1,12,14,74,180,2309)]
tableone_sa4 <- merge(tableone_sa4, PCE, by = "id", all.x = TRUE)
tableone_sa4 <- merge(tableone_sa4, LE8score_composite_, by = "id", all.x = TRUE)
tableone_sa4 <- merge(tableone_sa4, family, by = "id", all.x = TRUE)
tableone_sa4 <- merge(tableone_sa4, fh, by = "id", all.x = TRUE)

#remove repeated columns
tableone_sa4$Sex.y<-NULL
tableone_sa4$age.y<-NULL

colnames(tableone_sa4) <- c("id", "Sex", "Age", "White British Ancestry", "C Reactive Protein","Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                            "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")

tab7 <- CreateTableOne(vars = myVars, data = tableone_sa4, factorVars = catVars)
table1_sa4<-print(tab7,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_sa4, file = "table1_sa4.csv")