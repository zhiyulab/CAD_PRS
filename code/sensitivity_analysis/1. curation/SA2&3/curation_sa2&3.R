
library(data.table)
library(readr)

#######################
pheno <- fread(
  "/medpop/esp2/projects/UK_Biobank/baskets/2008463/ukb47823.csv.gz",
  select = c("eid", "20110-0.0", "20110-0.1","20110-0.2","20110-0.3","20110-0.4","20110-0.5","20110-0.6","20110-0.7","20110-0.8","20110-0.9","20110-0.10")
)

fwrite(pheno,"pheno_mom.csv")

cols <- paste0("20110-0.", 0:10)

ids_with_1 <- pheno_mom$eid[
  apply(pheno_mom[, cols], 1, function(x) any(x == 1, na.rm = TRUE))
]

pheno_dad <- read_csv("pheno_dad.csv")

cols <- paste0("20107-0.", 0:9)

ids_with_1_dad <- pheno_dad$eid[
  apply(pheno_dad[, cols], 1, function(x) any(x == 1, na.rm = TRUE))
]
family_history<-pheno_dad[,1]
family_history$mom<-ifelse(family_history$eid %in% ids_with_1,1,0)
family_history$dad<-ifelse(family_history$eid %in% ids_with_1_dad,1,0)
family_history$family<-ifelse(family_history$mom==1 | family_history$dad==1,1,0)
family_history$both<-ifelse(family_history$mom==1 & family_history$dad==1,1,0)
table(family_history$both)
family_history_both<-family_history[family_history$both==1,]
fwrite(family_history_both, "family_history_both.csv")

family_history_final$family_history<-ifelse(family_history_final$both==1,"both_parents",
                                            ifelse (family_history_final$mom==1,"mother",
                                                    ifelse(family_history_final$dad==1,"father","no_history")))
family_history_final<-family_history_final[,c(1,6)]
names(family_history_final)[1] <- "id"
fwrite(family_history_final, "family_history_final.csv")

#######################
# Load libraries
library(data.table)

# === Load second phenotype dataset for intermediate CAD ===
X20250615_pheno_noGP <- fread("/medpop/esp2/lli/UKB_pheno_curation/Results/20250615_pheno_noGP.tsv.gz")


# === Select individuals with intermediate CAD ===
ID_intermediate <- subset(X20250615_pheno_noGP, has_disease_Coronary_Artery_Disease_INTERMEDIATE == 1)

library(readr)
family_history<-read.csv("/home/unix/sliang/sliang/family/family_history_both.csv")
family_history<-family_history[family_history$id %in% ID_intermediate$sample_id,]
fwrite(family_history,"family_history_both_CAD.csv")