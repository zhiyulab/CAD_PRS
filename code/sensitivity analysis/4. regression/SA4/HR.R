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
outcome<-X20250615_pheno_noGP[,c(1,179,180,181)]
names(outcome)[1] <- "id"  
regression<-merge(regression,outcome,by = "id", all.x = TRUE)
regression<-regression[(regression$prevalent_disease_Coronary_Artery_Disease_INTERMEDIATE==0),]
regression<-regression[(!is.na(regression$incident_disease_Coronary_Artery_Disease_INTERMEDIATE)),]
regression<-regression[regression$followup_Coronary_Artery_Disease_INTERMEDIATE>0,]


###############################
#HR, with censoring
censor<-X20250615_pheno_noGP[,c(1,271,272,273,215,216,217,267,268,269)]
names(censor)[1] <- "id"  
regression2<-merge(regression,censor,by="id",all.x=TRUE)
regression2<-regression2[(regression2$prevalent_disease_Hypertension==0),]
regression2<-regression2[(regression2$prevalent_disease_Hypercholesterolemia==0),]
regression2<-regression2[(regression2$prevalent_disease_Diabetes_All==0),]
regression2<-regression2[regression2$followup_Diabetes_All>0,]
regression2<-regression2[regression2$followup_Hypercholesterolemia>0,]
regression2<-regression2[regression2$followup_Hypertension>0,]

#censor at the time of developing SMuRF (DM, hypertension or hypercholesterolmemia)
regression2$censor <- ifelse(
  (regression2$followup_Diabetes_All <= regression2$followup_Coronary_Artery_Disease_INTERMEDIATE |
     regression2$followup_Hypercholesterolemia <= regression2$followup_Coronary_Artery_Disease_INTERMEDIATE |
     regression2$followup_Hypertension <= regression2$followup_Coronary_Artery_Disease_INTERMEDIATE),
  1,
  ifelse(
    (regression2$incident_disease_Coronary_Artery_Disease_INTERMEDIATE == 0 &
       regression2$incident_disease_Diabetes_All == 0 &
       regression2$incident_disease_Hypercholesterolemia == 0 &
       regression2$incident_disease_Hypertension == 0),
    2,
    0
  )
)

regression2$followup <- ifelse(
  regression2$censor == 1,
  pmin(regression2$followup_Diabetes_All,
       regression2$followup_Hypercholesterolemia,
       regression2$followup_Hypertension), 
  ifelse(
    regression2$censor == 2,
    regression2$followup_Coronary_Artery_Disease_INTERMEDIATE,
    ifelse(
      regression2$censor == 0,
      regression2$followup_Coronary_Artery_Disease_INTERMEDIATE,
      NA  
    )
  )
)

regression2$event<-ifelse (regression2$censor!=0,0,1)
regression2_other<-regression2[regression2$in_white_British_ancestry_subset == 0, ]
regression2_white<-regression2[regression2$in_white_British_ancestry_subset == 1, ]


#all ancetry
library(survival)
library(tibble)
library(dplyr)
library(furrr)
library(future)

#Define exposure variables
exposures <- paste0("df", 1:58)

#Covariates
covariates <- "age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"  # edit as needed

run_cox <- function(var) {
  formula <- as.formula(
    paste0("Surv(followup, event) ~ ", var, " + ", covariates)
  )
  
  model <- tryCatch(
    coxph(formula, data = regression2),
    error = function(e) return(NULL)
  )
  
  if (is.null(model)) {
    return(tibble(Exposure = var, Term = NA, HR = NA, CI_lower = NA, CI_upper = NA, P_value = NA))
  }
  
  coefs <- summary(model)$coefficients
  rowname <- rownames(coefs)[1]
  
  hr <- exp(coefs[1, "coef"])
  se <- coefs[1, "se(coef)"]
  ci_low <- exp(coefs[1, "coef"] - 1.96 * se)
  ci_high <- exp(coefs[1, "coef"] + 1.96 * se)
  p <- coefs[1, "Pr(>|z|)"]
  
  tibble(
    Exposure = var,
    Term = rowname,
    HR = hr,
    CI_lower = ci_low,
    CI_upper = ci_high,
    P_value = p
  )
}

#Run all models in parallel
cox2_results <- future_map_dfr(exposures, run_cox, .progress = TRUE)

#white ancetry
#Define exposure variables
exposures <- paste0("df", 1:58)

#Covariates
covariates <- "age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"  # edit as needed

run_cox <- function(var) {
  formula <- as.formula(
    paste0("Surv(followup, event) ~ ", var, " + ", covariates)
  )
  
  model <- tryCatch(
    coxph(formula, data = regression2_white),
    error = function(e) return(NULL)
  )
  
  if (is.null(model)) {
    return(tibble(Exposure = var, Term = NA, HR = NA, CI_lower = NA, CI_upper = NA, P_value = NA))
  }
  
  coefs <- summary(model)$coefficients
  rowname <- rownames(coefs)[1]
  
  hr <- exp(coefs[1, "coef"])
  se <- coefs[1, "se(coef)"]
  ci_low <- exp(coefs[1, "coef"] - 1.96 * se)
  ci_high <- exp(coefs[1, "coef"] + 1.96 * se)
  p <- coefs[1, "Pr(>|z|)"]
  
  tibble(
    Exposure = var,
    Term = rowname,
    HR = hr,
    CI_lower = ci_low,
    CI_upper = ci_high,
    P_value = p
  )
}

#Run all models in parallel
cox2_results_white <- future_map_dfr(exposures, run_cox, .progress = TRUE)

#other ancestry
#Define exposure variables
exposures <- paste0("df", 1:58)

#Covariates
covariates <- "age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"  # edit as needed

run_cox <- function(var) {
  formula <- as.formula(
    paste0("Surv(followup, event) ~ ", var, " + ", covariates)
  )
  
  model <- tryCatch(
    coxph(formula, data = regression2_other),
    error = function(e) return(NULL)
  )
  
  if (is.null(model)) {
    return(tibble(Exposure = var, Term = NA, HR = NA, CI_lower = NA, CI_upper = NA, P_value = NA))
  }
  
  coefs <- summary(model)$coefficients
  rowname <- rownames(coefs)[1]
  
  hr <- exp(coefs[1, "coef"])
  se <- coefs[1, "se(coef)"]
  ci_low <- exp(coefs[1, "coef"] - 1.96 * se)
  ci_high <- exp(coefs[1, "coef"] + 1.96 * se)
  p <- coefs[1, "Pr(>|z|)"]
  
  tibble(
    Exposure = var,
    Term = rowname,
    HR = hr,
    CI_lower = ci_low,
    CI_upper = ci_high,
    P_value = p
  )
}

#Run all models in parallel
cox2_results_other <- future_map_dfr(exposures, run_cox, .progress = TRUE)
output_dir <- "/medpop/esp2/sliang"
fwrite(cox2_results,file.path(output_dir, "cox2_results.csv"))
fwrite(cox2_results_other,file.path(output_dir, "cox2_results_other.csv"))
fwrite(cox2_results_white,file.path(output_dir, "cox2_results_white.csv"))

