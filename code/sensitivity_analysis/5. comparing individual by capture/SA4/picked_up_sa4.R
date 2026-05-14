if (!require("tableone")) install.packages("tableone", repos="https://cloud.r-project.org")
library(tableone)
library(readr)

SMuRFless_55_sum <- read_csv("/medpop/esp2/sliang/SMuRFless_55_sum.csv")
LE8score_composite_ <- read_csv("/medpop/esp2/yxliu/LE8score_calculation/updated_LE8score_wComposite.txt")
pce_zy_2 <- read_delim("/medpop/esp2/sliang/pce_zy_2.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
SMuRFless_55 <- read_csv("/medpop/esp2/sliang/SMuRFless_CAD_under55.csv")
family<-read_csv("/home/unix/sliang/sliang/family/family_history_final.csv")
fh<-read_delim("/medpop/esp2/sliang/ukb_exome_450k_fh.carrier.txt", 
               delim = "\t", escape_double = FALSE, 
               trim_ws = TRUE)
fh<-fh[,c(1,3)]
colnames(fh)[1] <- "id"

#only include those with scores
SMuRFless_55_sum<-SMuRFless_55_sum[!is.na(SMuRFless_55_sum$df1_cutoff),]

#get columns of interest to present in descriptive table
SMuRFless_55_sum<-SMuRFless_55_sum[, c(1,60)]
onsetage<-SMuRFless_55[,c(1,2329)]
SMuRFless_55<-SMuRFless_55[,c(1,19:2328)]
ID<-merge(SMuRFless_55_sum,SMuRFless_55,by="id",all.x=TRUE)
ID<-merge(ID,LE8score_composite_,by="id",all.x=TRUE)
ID<-merge(ID,pce_zy_2,by="id",all.x=TRUE)
peopleprofile <- ID[, c(1,2,13,15,2310,2315,2316,2317,2318,2319,2320,2321,2322,2323,1392,21,647,939,2334,174,189,181)]
peopleprofile<-merge(peopleprofile,onsetage,by="id",all.x=TRUE)
peopleprofile<-merge(peopleprofile,family,by="id",all.x=TRUE)
peopleprofile<-merge(peopleprofile,fh,by="id",all.x=TRUE)
colnames(peopleprofile)

#rank by decreasing number of pickups
peopleprofile <- peopleprofile[order(-peopleprofile$sum), ]  
rownames(peopleprofile) <- 1:157# Descending

#those not picked up (cutoff=1) for any scores has sum=0
not_picked_up<-peopleprofile[(peopleprofile$sum==0),]
rownames(not_picked_up) <- 1:26

#get a similar size population of frequently picked up individuals 
picked_up<-peopleprofile[c(1:22),]

colnames(picked_up) <- c("id", "sum", "Sex", "Age", "Townsend Deprivation Index", "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                         "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score", "Race", "Alcohol intake frequency", "Chronic Kidney Disease",
                         "Rheumatoid arthritis", "Pooled Cohort Equation Risk","Apolipoprotein B", "Lipoprotein A", "C Reactive Protein", "Age of CAD Onset","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
colnames(not_picked_up) <- c("id", "sum", "Sex", "Age", "Townsend Deprivation Index", "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
                             "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score", "Race", "Alcohol intake frequency", "Chronic Kidney Disease",
                             "Rheumatoid arthritis", "Pooled Cohort Equation Risk","Apolipoprotein B", "Lipoprotein A", "C Reactive Protein", "Age of CAD Onset","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
write.csv(picked_up, file = "picked_up_sa4.csv", row.names = FALSE)
write.csv(not_picked_up, file = "not_picked_up_sa4.csv", row.names = FALSE)


## Vector of variables to summarize
myVars <- c("Sex", "Age", "Townsend Deprivation Index", "LE8 Diet Score", "LE8 PA Score", "LE8 Smoking Score", "LE8 Sleep Health Points", "LE8 BMI Points",
            "LE8 Blood Lipid Points", "LE8 HbA1c Score", "LE8 BP Score", "LE8 Composite Score", "Race", "Alcohol intake frequency", "Chronic Kidney Disease",
            "Rheumatoid arthritis", "Pooled Cohort Equation Risk","Apolipoprotein B","Lipoprotein A", "Age of CAD Onset","C Reactive Protein","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
## Vector of categorical variables that need transformation
catVars <- c("Sex", "Race", "Alcohol intake frequency", "Chronic Kidney Disease",
             "Rheumatoid arthritis","Family Heart Disease History","Familial Hypercholesterolemia Variant Status")
## Create a TableOne object
tab1 <- CreateTableOne(vars = myVars, data = not_picked_up, factorVars = catVars)
table1_notpickedup<-print(tab1,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_notpickedup, file = "Table1_notpickedup.csv")

tab2 <- CreateTableOne(vars = myVars, data = picked_up, factorVars = catVars)
table1_pickedup<-print(tab2,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1_pickedup, file = "Table1_pickedup.csv")

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
write.csv(pval_continuous, file = "pval_continuous.csv", row.names = FALSE)

#get p value for categorical variables with chi-square or fisher exact
categorical_cols <-c("Sex", "Race", "Alcohol intake frequency", "Chronic Kidney Disease",
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

write.csv(pval_categorical, file = "pval_categorical.csv")