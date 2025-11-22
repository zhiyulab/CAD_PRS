if (!require("tableone")) install.packages("tableone", repos="https://cloud.r-project.org")
library(tableone)
if (!require("readr")) install.packages("readr", repos = "https://cloud.r-project.org")
library(readr)
LE8score_composite_ <- read_csv("/medpop/esp2/yxliu/LE8score_calculation/updated_LE8score_wComposite.txt")
pce_zy_2 <- read_delim("/medpop/esp2/sliang/pce_zy_2.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
PCE <- pce_zy_2[, c(1,12)]

#I put these files (output from step 1) in the code check folder but they are also available in the following directory
SMuRF<- read_csv("/medpop/esp2/sliang/SMuRF.csv")
SMuRFless_55 <- read_csv("/medpop/esp2/sliang/SMuRFless_CAD_under55.csv")
SMuRFless_all <- read_csv("/medpop/esp2/sliang/SMuRFless_CAD.csv")

#get columns with variables of interest to present in table one
tableone55<-SMuRFless_55[,c(1, 19:2329)]
tableone55<-tableone55[, c(1,12,14,74,180,2309)]
tableone55 <- merge(tableone55, PCE, by = "id", all.x = TRUE)
tableone55 <- merge(tableone55, LE8score_composite_, by = "id", all.x = TRUE)
#remove repeated columns
tableone55$Sex.y<-NULL
tableone55$age.y<-NULL
colnames(tableone55) <- c("id", "Sex", "Age", "White British Ancestry", "C Reactive Protein","Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                             "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                             "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score")

#repeat for SMuRFless group
tableoneSMuRFless<-SMuRFless_all[,c(1, 19:2329)]
tableoneSMuRFless<-tableoneSMuRFless[, c(1,12,14,74,180,2309)]
tableoneSMuRFless <- merge(tableoneSMuRFless, PCE, by = "id", all.x = TRUE)
tableoneSMuRFless <- merge(tableoneSMuRFless, LE8score_composite_, by = "id", all.x = TRUE)
tableoneSMuRFless$Sex.y<-NULL
tableoneSMuRFless$age.y<-NULL
colnames(tableoneSMuRFless) <- c("id", "Sex", "Age", "White British Ancestry","C Reactive Protein", "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                          "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                          "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score")

#repeat for SMuRF group
tableoneSMuRF <- SMuRF[, c(1,12,14,74,180,2309)]
tableoneSMuRF <- merge(tableoneSMuRF, PCE, by = "id", all.x = TRUE)
tableoneSMuRF <- merge(tableoneSMuRF, LE8score_composite_, by = "id", all.x = TRUE)
tableoneSMuRF$Sex.y<-NULL
tableoneSMuRF$age.y<-NULL
colnames(tableoneSMuRF) <- c("id", "Sex", "Age", "White British Ancestry","C Reactive Protein", "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                             "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                             "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score")



## Vector of variables to summarize
myVars <- c("Sex", "Age", "White British Ancestry", "C Reactive Protein", "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
            "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score")
## Vector of categorical variables that need transformation
catVars <- c("Sex", "White British Ancestry")
## Create a TableOne object
tab1 <- CreateTableOne(vars = myVars, data = tableone55, factorVars = catVars)
table1_55<-print(tab1,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_55, file = "Table1_55.csv")

tab2 <- CreateTableOne(vars = myVars, data = tableoneSMuRFless, factorVars = catVars)
table1_SMuRFless<-print(tab2,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_SMuRFless, file = "Table1_SMuRFless.csv")

tab3 <- CreateTableOne(vars = myVars, data = tableoneSMuRF, factorVars = catVars)
table1_SMuRF<-print(tab3,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_SMuRF, file = "Table1_SMuRF.csv")

#get p value for continuous variables with t test
continuous_cols <- c("Age","C Reactive Protein", "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                     "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                     "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score")  # replace with your variables
p_values <- sapply(continuous_cols, function(col) {
  t.test(tableoneSMuRFless[[col]], tableoneSMuRF[[col]])$p.value
})

pval_continuous_SMuRF <- data.frame(
  variable = continuous_cols,
  p_value = as.numeric(p_values)  # Remove names to avoid confusion
)
write.csv(pval_continuous_SMuRF, file = "pval_continuous_SMuRF.csv", row.names = FALSE)

#get p value for continuous variables with t test
continuous_cols <- c("Age","C Reactive Protein", "Townsend Deprivation Index", "Pooled Cohort Equation Risk",
                     "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                     "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score")  # replace with your variables
p_values <- sapply(continuous_cols, function(col) {
  t.test(tableoneSMuRFless[[col]], tableone55[[col]])$p.value
})

pval_continuous_age <- data.frame(
  variable = continuous_cols,
  p_value = as.numeric(p_values)  # Remove names to avoid confusion
)
write.csv(pval_continuous_age, file = "pval_continuous_age.csv", row.names = FALSE)

#get p value for categorical variables with chi-square or fisher exact
categorical_cols <-c("Sex","White British Ancestry")
get_pvals <- function(col) {
  df <- rbind(
    data.frame(Group = "SMuRFless", Var = tableoneSMuRFless[[col]]),
    data.frame(Group = "SMuRF", Var = tableoneSMuRF[[col]])
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
pval_categorical_SMuRF <- data.frame(
  variable = rownames(pval_results),
  chi_sq_p = pval_results[, "chi_sq_p"],
  fisher_p = pval_results[, "fisher_p"],
  row.names = NULL
)

write.csv(pval_categorical_SMuRF, file = "pval_categorical_SMuRF.csv")

#get p value for categorical variables with chi-square or fisher exact
categorical_cols <-c("Sex","White British Ancestry")
get_pvals <- function(col) {
  df <- rbind(
    data.frame(Group = "SMuRFless", Var = tableoneSMuRFless[[col]]),
    data.frame(Group = "earlyonset", Var = tableone55[[col]])
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
pval_categorical_55 <- data.frame(
  variable = rownames(pval_results),
  chi_sq_p = pval_results[, "chi_sq_p"],
  fisher_p = pval_results[, "fisher_p"],
  row.names = NULL
)

write.csv(pval_categorical_55, file = "pval_categorical_55.csv")