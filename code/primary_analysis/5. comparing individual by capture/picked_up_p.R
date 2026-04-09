if (!require("tableone")) install.packages("tableone", repos="https://cloud.r-project.org")
library(tableone)
library(readr)
library(data.table)

cohort_p_sum <- read.csv("/medpop/esp2/sliang/cohort_p_sum.csv", na.strings = c("", " "))

cohort_p <- read_csv("/medpop/esp2/sliang/cohort_p_55.csv")
LE8score_composite_ <- read_csv("/medpop/esp2/yxliu/LE8score_calculation/updated_LE8score_wComposite.txt")
pce_zy_2 <- read_delim("/medpop/esp2/sliang/pce_zy_2.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
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

#only include those with scores

cohort_p_sum<-cohort_p_sum[!is.na(cohort_p_sum$df1_cutoff),]

#get columns of interest to present in descriptive table
cohort_p_sum<-cohort_p_sum[, c(1,60)]
onsetage<-cohort_p[,c(1,20)]

ID<-merge(cohort_p_sum,pheno,by="id",all.x=TRUE)
ID<-merge(ID,LE8score_composite_,by="id",all.x=TRUE)
ID<-merge(ID,pce_zy_2,by="id",all.x=TRUE)
peopleprofile <- ID[, c(1,2,13,15,16,2310,2315,2316,2317,2318,2319,2320,2321,2322,2323,1392,21,647,939,2334,174,189,181)]
peopleprofile<-merge(peopleprofile,onsetage,by="id",all.x=TRUE)
peopleprofile<-merge(peopleprofile,family,by="id",all.x=TRUE)
peopleprofile<-merge(peopleprofile,fh,by="id",all.x=TRUE)

colnames(peopleprofile)

#rank by decreasing number of pickups
peopleprofile <- peopleprofile[order(-peopleprofile$sum), ]  
rownames(peopleprofile) <- 1:nrow(peopleprofile)# Descending

#those not picked up (cutoff=1) for any scores has sum=0
not_picked_up<-peopleprofile[(peopleprofile$sum==0),]
nrow(not_picked_up)
rownames(not_picked_up) <- NULL
rownames(not_picked_up) <- 1:nrow(not_picked_up)

#get a similar size population of frequently picked up individuals 
picked_up<-peopleprofile[(peopleprofile$sum>20),]

colnames(picked_up) <- c("id", "sum", "Sex", "Age", "Smoking","Townsend Deprivation Index", "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                         "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score", "Race", "Alcohol intake frequency", "Chronic Kidney Disease",
                         "Rheumatoid arthritis", "Pooled Cohort Equation Risk","Apolipoprotein B", "Lipoprotein A", "C Reactive Protein", "Age of CAD Onset","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
colnames(not_picked_up) <- c("id", "sum", "Sex", "Age", "Smoking", "Townsend Deprivation Index", "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                             "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score", "Race", "Alcohol intake frequency", "Chronic Kidney Disease",
                             "Rheumatoid arthritis", "Pooled Cohort Equation Risk","Apolipoprotein B", "Lipoprotein A", "C Reactive Protein", "Age of CAD Onset","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")

## Vector of variables to summarize
myVars <- c("Sex", "Age", "Smoking","Townsend Deprivation Index", "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score", "Race", "Alcohol intake frequency", "Chronic Kidney Disease",
            "Rheumatoid arthritis", "Pooled Cohort Equation Risk","Apolipoprotein B","Lipoprotein A", "Age of CAD Onset","C Reactive Protein","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
## Vector of categorical variables that need transformation
catVars <- c("Sex", "Smoking","Race", "Alcohol intake frequency", "Chronic Kidney Disease",
             "Rheumatoid arthritis","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
## Create a TableOne object
tab1 <- CreateTableOne(vars = myVars, data = not_picked_up, factorVars = catVars)
table1_notpickedup<-print(tab1,quote = FALSE, noSpaces = TRUE, printToggle = FALSE,showAllLevels = TRUE, pDigits = 3)
write.csv(table1_notpickedup, file = "Table1_notpickedup_p.csv")

tab2 <- CreateTableOne(vars = myVars, data = picked_up, factorVars = catVars)
table1_pickedup<-print(tab2,quote = FALSE, noSpaces = TRUE, printToggle = FALSE,showAllLevels = TRUE, pDigits = 3)
write.csv(table1_pickedup, file = "Table1_pickedup_p.csv")

#get p value for continuous variables with t test
continuous_cols <- c("Age", "Townsend Deprivation Index", "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                     "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score", "Pooled Cohort Equation Risk","Apolipoprotein B","Lipoprotein A","C Reactive Protein","Age of CAD Onset")  # replace with your variables
p_values <- sapply(continuous_cols, function(col) {
  t.test(not_picked_up[[col]], picked_up[[col]])$p.value
})

pval_continuous <- data.frame(
  variable = continuous_cols,
  p_value = as.numeric(p_values)  # Remove names to avoid confusion
)
write.csv(pval_continuous, file = "pval_continuous_p.csv", row.names = FALSE)

#get p value for categorical variables with chi-square or fisher exact
categorical_cols <-c("Sex", "Smoking","Race", "Alcohol intake frequency", "Chronic Kidney Disease",
                     "Rheumatoid arthritis","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
get_pvals <- function(col) {
  df <- rbind(
    data.frame(Group = "Not_Picked_Up", Var = not_picked_up[[col]]),
    data.frame(Group = "Picked_Up", Var = picked_up[[col]])
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
pval_categorical <- data.frame(
  variable = rownames(pval_results),
  chi_sq_p = pval_results[, "chi_sq_p"],
  fisher_p = pval_results[, "fisher_p"],
  row.names = NULL
)

write.csv(pval_categorical, file = "pval_categorical_p.csv")
