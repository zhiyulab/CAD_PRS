
library(data.table)
library(readr)
pheno<-fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt") 
pheno = data.frame(pheno)
pheno$IDs_toRemove_SampleQCv3 = ifelse((!(pheno$Submitted_Gender == pheno$Inferred_Gender) | pheno$Non_Consented== 1),1,0) 
pheno = pheno[which(pheno$IDs_toRemove_SampleQCv3==0),]
UKB_PRS<- read_csv("/home/unix/sliang/sliang/UKB_PRS.csv")
X20250615_pheno_noGP <- fread("/medpop/esp2/lli/UKB_pheno_curation/Results/20250615_pheno_noGP.tsv.gz")

pheno1<-pheno[, c(1:14)]
names(UKB_PRS)[1] <- "id"  
regression<-merge(pheno1,UKB_PRS,by = "id", all.y = TRUE)
race<-pheno[,c("id","in_white_British_ancestry_subset")]
regression<-merge(regression,race,by = "id", all.x = TRUE)
outcome<-X20250615_pheno_noGP[,c(1,178,179,180,181)]
names(outcome)[1] <- "id"  
regression<-merge(regression,outcome,by = "id", all.x = TRUE)
write.table(regression, "regression.txt", sep = "\t", row.names = FALSE)

regression_other<-regression[regression$in_white_British_ancestry_subset == 0, ]
regression_white<-regression[regression$in_white_British_ancestry_subset == 1, ]

###############################
#OR

library(speedglm)
library(purrr)
library(tibble)

# Define your exposure variables
exposures <- paste0("df", 1:58)

#all ancestry
run_logistic_fast <- function(var) {
  formula <- as.formula(paste("has_disease_Coronary_Artery_Disease_INTERMEDIATE ~", var, "+ age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))  # include other covariates
  model <- speedglm(formula, data = regression, family = binomial())
  coefs <- summary(model)$coefficients
  
  # Extract the row for the exposure
  or <- exp(coefs[var, "Estimate"])
  se <- coefs[var, "Std. Error"]
  ci_low <- exp(coefs[var, "Estimate"] - 1.96 * se)
  ci_high <- exp(coefs[var, "Estimate"] + 1.96 * se)
  p <- coefs[var, "Pr(>|z|)"]
  
  tibble(
    Exposure = var,
    OR = or,
    CI_lower = ci_low,
    CI_upper = ci_high,
    P_value = p
  )
}

results <- map_dfr(exposures, run_logistic_fast)

#white ancestry
run_logistic_fast <- function(var) {
  formula <- as.formula(paste("has_disease_Coronary_Artery_Disease_INTERMEDIATE ~", var, "+ age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))  # include other covariates
  model <- speedglm(formula, data = regression_white, family = binomial())
  coefs <- summary(model)$coefficients
  
  # Extract the row for the exposure
  or <- exp(coefs[var, "Estimate"])
  se <- coefs[var, "Std. Error"]
  ci_low <- exp(coefs[var, "Estimate"] - 1.96 * se)
  ci_high <- exp(coefs[var, "Estimate"] + 1.96 * se)
  p <- coefs[var, "Pr(>|z|)"]
  
  tibble(
    Exposure = var,
    OR = or,
    CI_lower = ci_low,
    CI_upper = ci_high,
    P_value = p
  )
}

results_white <- map_dfr(exposures, run_logistic_fast)

#other ancestry
run_logistic_fast <- function(var) {
  formula <- as.formula(paste("has_disease_Coronary_Artery_Disease_INTERMEDIATE ~", var, "+ age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))  # include other covariates
  model <- speedglm(formula, data = regression_other, family = binomial())
  coefs <- summary(model)$coefficients
  
  # Extract the row for the exposure
  or <- exp(coefs[var, "Estimate"])
  se <- coefs[var, "Std. Error"]
  ci_low <- exp(coefs[var, "Estimate"] - 1.96 * se)
  ci_high <- exp(coefs[var, "Estimate"] + 1.96 * se)
  p <- coefs[var, "Pr(>|z|)"]
  
  tibble(
    Exposure = var,
    OR = or,
    CI_lower = ci_low,
    CI_upper = ci_high,
    P_value = p
  )
}

results_other <- map_dfr(exposures, run_logistic_fast)

output_dir <- "/medpop/esp2/sliang"
fwrite(results,file.path(output_dir, "results.csv"))
fwrite(results_other,file.path(output_dir, "results_other.csv"))
fwrite(results_white,file.path(output_dir, "results_white.csv"))

