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

#I put these files (output from step 1) in the code check folder but they are also available in the following directory
cohort_p_all<- read_csv("/medpop/esp2/sliang/cohort_p_all.csv")
cohort_p_55<- read_csv("/medpop/esp2/sliang/cohort_p_55.csv")
not_cohort_p<- read_csv("/medpop/esp2/sliang/not_cohort_p.csv")
cohort_p_all<-cohort_p_all[,"id"]
cohort_p_55<-cohort_p_55[,"id"]
not_cohort_p<-not_cohort_p[,"id"]

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

#get columns with variables of interest to present in table one
tableone_55<-merge(cohort_p_55, pheno, by = "id", all.x = TRUE)
tableone_55<-tableone_55[, c(1,12,14,74,180,2309)]
tableone_55 <- merge(tableone_55, PCE, by = "id", all.x = TRUE)
tableone_55 <- merge(tableone_55, LE8score_composite_, by = "id", all.x = TRUE)
tableone_55 <- merge(tableone_55, family, by = "id", all.x = TRUE)
tableone_55 <- merge(tableone_55, fh, by = "id", all.x = TRUE)

#remove repeated columns
tableone_55$Sex.y<-NULL
tableone_55$age.y<-NULL
colnames(tableone_55)

#get columns with variables of interest to present in table one
tableone_all<-merge(cohort_p_all, pheno, by = "id", all.x = TRUE)
tableone_all<-tableone_all[, c(1,12,14,74,180,2309)]
tableone_all <- merge(tableone_all, PCE, by = "id", all.x = TRUE)
tableone_all <- merge(tableone_all, LE8score_composite_, by = "id", all.x = TRUE)
tableone_all <- merge(tableone_all, family, by = "id", all.x = TRUE)
tableone_all <- merge(tableone_all, fh, by = "id", all.x = TRUE)

#remove repeated columns
tableone_all$Sex.y<-NULL
tableone_all$age.y<-NULL

#get columns with variables of interest to present in table one
tableone_not<-merge(not_cohort_p, pheno, by = "id", all.x = TRUE)
tableone_not<-tableone_not[, c(1,12,14,74,180,2309)]
tableone_not <- merge(tableone_not, PCE, by = "id", all.x = TRUE)
tableone_not <- merge(tableone_not, LE8score_composite_, by = "id", all.x = TRUE)
tableone_not <- merge(tableone_not, family, by = "id", all.x = TRUE)
tableone_not <- merge(tableone_not, fh, by = "id", all.x = TRUE)

#remove repeated columns
tableone_not$Sex.y<-NULL
tableone_not$age.y<-NULL

colnames(tableone_not) <- c("id", "Sex", "Age", "White British Ancestry", "C Reactive Protein","Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                             "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                             "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
colnames(tableone_all) <- c("id", "Sex", "Age", "White British Ancestry", "C Reactive Protein","Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                            "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
colnames(tableone_55) <- c("id", "Sex", "Age", "White British Ancestry", "C Reactive Protein","Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                            "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")

## Vector of variables to summarize
myVars <- c("Sex", "Age", "White British Ancestry", "C Reactive Protein", "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
            "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
## Vector of categorical variables that need transformation
catVars <- c("Sex", "White British Ancestry","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
## Create a TableOne object
tab1 <- CreateTableOne(vars = myVars, data = tableone_all, factorVars = catVars)
table1_pall<-print(tab1,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_pall, file = "Table1_pall.csv")

## Create a TableOne object
tab2 <- CreateTableOne(vars = myVars, data = tableone_55, factorVars = catVars)
table1_p55<-print(tab2,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_p55, file = "Table1_p55.csv")

## Create a TableOne object
tab3 <- CreateTableOne(vars = myVars, data = tableone_not, factorVars = catVars)
table1_pnot<-print(tab3,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_pnot, file = "table1_pnot.csv")

#get p value for continuous variables with t test
continuous_cols <- c("Age","C Reactive Protein", "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                     "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                     "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score")  # replace with your variables
# Get raw p-values
p_values <- sapply(continuous_cols, function(col) {
  t.test(tableone_all[[col]], tableone_not[[col]])$p.value
})

# Format extremely small values instead of showing 0
p_values_formatted <- format.pval(p_values, eps = 1e-300)

pval_continuous_risk <- data.frame(
  variable = continuous_cols,
  p_value = p_values_formatted
)

write.csv(pval_continuous_risk, file = "pval_continuous_risk.csv", row.names = FALSE)

#get p value for continuous variables with t test
continuous_cols <- c("Age","C Reactive Protein", "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                     "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                     "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score")  # replace with your variables
p_values <- sapply(continuous_cols, function(col) {
  t.test(tableone_55[[col]], tableone_all[[col]])$p.value
})

pval_continuous_age <- data.frame(
  variable = continuous_cols,
  p_value = as.numeric(p_values)  # Remove names to avoid confusion
)
write.csv(pval_continuous_age, file = "pval_continuous_page.csv", row.names = FALSE)

#get p value for categorical variables with chi-square or fisher exact
categorical_cols <-c("Sex","White British Ancestry","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
get_pvals <- function(col) {
  df <- rbind(
    data.frame(Group = "highrisk", Var = tableone_all[[col]]),
    data.frame(Group = "nothighrisk", Var = tableone_not[[col]])
  )
  table_data <- table(df$Group, df$Var)
  # Try fisher.test with error handling
  fisher_p <- tryCatch(
    fisher.test(table_data)$p.value,
    error = function(e) NA
  )
  # Try chisq.test with error handling
  chisq_p <- tryCatch(
    chisq.test(table_data)$p.value,
    error = function(e) NA
  )
  return(c(chi_sq_p = chisq_p, fisher_p = fisher_p))
}
# Apply over all categorical variables
pval_results <- t(sapply(categorical_cols, get_pvals))
# Convert to dataframe
pval_categorical_risk <- data.frame(
  variable = rownames(pval_results),
  chi_sq_p = pval_results[, "chi_sq_p"],
  fisher_p = pval_results[, "fisher_p"],
  row.names = NULL
)

write.csv(pval_categorical_risk, file = "pval_categorical_risk.csv")

#get p value for categorical variables with chi-square or fisher exact
categorical_cols <-c("Sex","White British Ancestry","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
get_pvals <- function(col) {
  df <- rbind(
    data.frame(Group = "all", Var = tableone_all[[col]]),
    data.frame(Group = "earlyonset", Var = tableone_55[[col]])
  )
  table_data <- table(df$Group, df$Var)
  # Try fisher.test with error handling
  fisher_p <- tryCatch(
    fisher.test(table_data)$p.value,
    error = function(e) NA
  )
  # Try chisq.test with error handling
  chisq_p <- tryCatch(
    chisq.test(table_data)$p.value,
    error = function(e) NA
  )
  return(c(chi_sq_p = chisq_p, fisher_p = fisher_p))
}
# Apply over all categorical variables
pval_results <- t(sapply(categorical_cols, get_pvals))
# Convert to dataframe
pval_categorical_age <- data.frame(
  variable = rownames(pval_results),
  chi_sq_p = pval_results[, "chi_sq_p"],
  fisher_p = pval_results[, "fisher_p"],
  row.names = NULL
)

write.csv(pval_categorical_age, file = "pval_categorical_page.csv")
