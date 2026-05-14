library(data.table)
library(readr)



####################################################################################

MGBB<-MGBB_Phenos_2CODE_2022_06_21[MGBB_Phenos_2CODE_2022_06_21$Has_cad==1,]
MGBB<-MGBB[,c(1,2,3,4,5,8,9,10,11,12,15,17,20,27,30,37,40)]

library(dplyr)
library(lubridate)

MGBB_filtered <- MGBB %>%
  filter(
    is.na(dm_Date) | 
      is.na(cad_Date) | 
      dm_Date >= cad_Date %m+% months(6)
  )


MGBB_filtered <- MGBB_filtered %>%
  filter(
    is.na(htn_Date) | 
      is.na(cad_Date) | 
      htn_Date >= cad_Date %m+% months(6)
  )

MGBB_filtered <- MGBB_filtered %>%
  filter(
    is.na(hyperlipidemia_Date) | 
      is.na(cad_Date) | 
      hyperlipidemia_Date >= cad_Date %m+% months(6)
  )

IDs<-MGBB_filtered[,c(2,3,5,14)]

#############################################################
IDs_postphy<-IDs[IDs$phy==0,]

lab_6<-fread("PN16_20250218_110101-1_Lab.txt")
lab_6<- lab_6[!is.na(match(lab_6$EMPI, IDs_postphy$EMPI))]

include<- "(LDL|CHOL|GHBA1C)"
out<-"(CHOL/HDL|VLDL)"
lab_6 <- lab_6[
  grepl(include, Group_Id, ignore.case = TRUE) &
    !grepl(out, Group_Id, ignore.case = TRUE)
]
  
cad_dates<-IDs_postphy[,c(1,4)]
lab_6 <- merge(lab_6, cad_dates, by = "EMPI", all.x = TRUE)
setDT(lab_6)
lab_6[, Seq_Date_Time := as.Date(Seq_Date_Time, format = "%m/%d/%Y %H:%M")]
lab_6[, cad_Date := as.Date(cad_Date)]

# Compute distance
lab_6[, diff_days := as.numeric(difftime(Seq_Date_Time, cad_Date, units = "days"))]
lab_6$Group_Id[grepl("LDL", lab_6$Group_Id, ignore.case = TRUE)] <- "LDL"


lab_6 <- lab_6[
  ,
  {
    before <- .SD[diff_days < 0]
    
    if (nrow(before) > 0) {
      # choose closest BEFORE (largest negative)
      before[which.max(diff_days)]
    } else {
      after <- .SD[diff_days > 0]
      if (nrow(after) > 0) {
        # choose closest AFTER (smallest positive)
        after[which.min(diff_days)]
      } else {
        NULL
      }
    }
  },
  by = .(EMPI, Group_Id)
]


lab_6<-lab_6[,c(1,6,9,10,13,23,24)]

lab_all <- rbind(lab_1,lab_2,lab_3,lab_4,lab_5,lab_6)
rm(lab_1,lab_2,lab_3,lab_4,lab_5,lab_6)

lab_all$Test_Description[grepl("LDL", lab_all$Test_Description, ignore.case = TRUE)] <- "LDL"
lab_all$Test_Description[grepl("CHOL", lab_all$Test_Description, ignore.case = TRUE)] <- "CHOL"
lab_all[
  !(Test_Description %in% c("CHOL", "LDL")),
  Test_Description := "HBA1C"
]

hba1c <- lab_all[lab_all$Test_Description != "LDL" & lab_all$Test_Description != "CHOL",]
lipid<-lab_all[lab_all$Test_Description != "HBA1C",]

dm_date<-MGBB_Phenos_2CODE_2022_06_21[!is.na(match(MGBB_Phenos_2CODE_2022_06_21$EMPI, hba1c$EMPI)),]
lipid_date<-MGBB_Phenos_2CODE_2022_06_21[!is.na(match(MGBB_Phenos_2CODE_2022_06_21$EMPI, lipid$EMPI)),]
dm_date<-dm_date[,c(2,17)]
lipid_date<-lipid_date[,c(2,37)]
hba1c<-merge(hba1c,dm_date,by="EMPI")
lipid<-merge(lipid,lipid_date,by="EMPI")

hba1c_keep <- hba1c[
  is.na(hba1c$dm_Date) | hba1c$dm_Date > hba1c$Seq_Date_Time,
]
lipid_keep <- lipid[
  is.na(lipid$hyperlipidemia_Date) | lipid$hyperlipidemia_Date > lipid$Seq_Date_Time,
]
hba1c_keep$dm_Date<-NULL
lipid_keep$hyperlipidemia_Date<-NULL
lab_all_final<-rbind(hba1c_keep,lipid_keep)

setDT(lab_all_final)
setorder(lab_all_final,EMPI,Seq_Date_Time)
lab_all_final<-lab_all_final[,c(1,3,4)]
wide_lab <- lab_all_final %>%
  pivot_wider(
    id_cols = c(EMPI),
    names_from = Test_Description,
    values_from = Result
  )

wide_lab$LDL<-as.numeric(wide_lab$LDL)/38.67
wide_lab$CHOL<-as.numeric(wide_lab$CHOL)/38.67

wide_lab$IDs_toRemove <- ifelse(
  (!is.na(wide_lab$HBA1C) & wide_lab$HBA1C >= 6.5) |
    (!is.na(wide_lab$CHOL) & wide_lab$CHOL > 5.5) |
    (!is.na(wide_lab$LDL) & wide_lab$LDL > 3.5),
  1, 0
)

wide_lab_final<-wide_lab[,c(1,5)]
IDs<-merge(IDs,wide_lab_final,by="EMPI",all.x=T)
IDs$IDs_toRemove[is.na(IDs$IDs_toRemove)] <- 0
names(IDs)[names(IDs) == "IDs_toRemove"] <- "lab"


##########################################################
Phy_6<-fread("PN16_20250218_110101-6_Phy.txt")

Phy_6<- Phy_6[!is.na(match(Phy_6$EMPI, IDs$EMPI))]

include<- "(systolic|diastolic|cigarette|smok)"
out<-"(systolic/diastolic|pressure|ready|smokeless)"

# Filter
Phy_6 <- Phy_6[
  grepl(include, Concept_Name, ignore.case = TRUE) &
    !grepl(out, Concept_Name, ignore.case = TRUE)
]

Phy_6 <- merge(Phy_6, IDs, by = "EMPI", all.x = TRUE)
setDT(Phy_6)
Phy_6[, Date := as.Date(Date, format = "%m/%d/%Y")]
Phy_6[, cad_Date := as.Date(cad_Date)]

# Compute distance
Phy_6[, diff_days := as.numeric(difftime(Date, cad_Date, units = "days"))]

# Keep only closest test(s) for each empi and test group
Phy_6$Concept_Name[grepl("diastolic", Phy_6$Concept_Name, ignore.case = TRUE)] <- "Diastolic"
Phy_6$Concept_Name[grepl("systolic", Phy_6$Concept_Name, ignore.case = TRUE)] <- "Systolic"

Phy_6 <- Phy_6[
  ,
  {
    before <- .SD[diff_days < 0]
    after  <- .SD[diff_days > 0]
    
    rbind(
      before[diff_days == max(diff_days)],
      after[diff_days == min(diff_days)]
    )
  },
  by = .(EMPI, Concept_Name)
]


Phy_6<-Phy_6[,c(1,2,6,9,18,19)]
Phy_5<-Phy_6

phy_all <- rbind(Phy_1,Phy_2,Phy_3,Phy_4,Phy_5,Phy_6)
fwrite(phy_all,"phy_all.csv")
unique_phy_ID<-as.data.frame(unique(phy_all$EMPI))
fwrite(unique_phy_ID,"unique_phy_ID.csv")


include<- "(smok)"

# Filter
phy_all <- phy_all[
  grepl(include, phy_all$Concept_Name, ignore.case = TRUE),]

library(data.table)

# ensure data.table
setDT(phy_all)

phy_all[
  grepl("^Smoking Tobacco Use-", Concept_Name) & is.na(Result),
  `:=`(
    Result = sub("^Smoking Tobacco Use-", "", Concept_Name),
    Concept_Name = "Smoking Tobacco Use"
  )
]

phy_all[
  grepl("^Smoking (Start|Quit) Date", Concept_Name),
  `:=`(
    Date = as.IDate(Result),
    Result = fifelse(
      grepl("^Smoking Start Date", Concept_Name),
      "Start",
      "Quit"
    )
  )
]

phy_all[, diff_days := as.numeric(difftime(Date, cad_Date, units = "days"))]

# Keep only closest test(s) for each empi and test group
phy_all$Concept_Name[grepl("smoking", phy_all$Concept_Name, ignore.case = TRUE)] <- "smoking"


table(phy_all$Result)
setorder(phy_all, EMPI, Date)

phy_all$Result[grepl("unknown", phy_all$Result, ignore.case = TRUE)] <- "unknown"
phy_all$Result[grepl("yes|start|current", phy_all$Result, ignore.case = TRUE)] <- "smoking"
phy_all$Result[grepl("former", phy_all$Result, ignore.case = TRUE)] <- "former"
phy_all$Result[grepl("assessed", phy_all$Result, ignore.case = TRUE)] <- "unknown"
phy_all$Result[grepl("never|non", phy_all$Result, ignore.case = TRUE)] <- "non-smoker"

fwrite(phy_all,"phy_all_smoking.csv")

table(phy_all_smoking$Result)
phy_all_smoking$Result[grepl("quit", phy_all_smoking$Result, ignore.case = TRUE)] <- "quit"

setDT(phy_all_smoking)
phy_all_smoking[, diff_days := as.numeric(difftime(Date, cad_Date, units = "days"))]


phy_all_smoking <- phy_all_smoking[
  ,
  {
    before <- .SD[diff_days < 0]
    after  <- .SD[diff_days > 0]
    
    rbind(
      before[diff_days == max(diff_days)],
      after[diff_days == min(diff_days)]
    )
  },
  by = .(EMPI)
]
phy_all_smoking[diff_days > 0, Concept_Name := "smoking_after"]
phy_all_smoking[diff_days < 0, Concept_Name := "smoking_before"]


library(tidyr)
library(dplyr)

wide_phy_smok <- phy_all_smoking %>%
  pivot_wider(
    id_cols = EMPI,
    names_from = Concept_Name,
    values_from = Result
  )
wide_phy_smok$smoking_before<-as.character(wide_phy_smok$smoking_before)
wide_phy_smok$smoking_after<-as.character(wide_phy_smok$smoking_after)

fwrite(wide_phy_smok,"phy_all_smoking_final.csv")


smok<-phy_all_smoking_final[,c(1,4)]

htn_date<-MGBB_Phenos_2CODE_2022_06_21[!is.na(match(MGBB_Phenos_2CODE_2022_06_21$EMPI, phy_all$EMPI)),]
htn_date<-htn_date[,c(2,12)]
phy_all<-merge(phy_all,htn_date,by="EMPI")
htn <- phy_all[
  phy_all$Concept_Name == "Diastolic" | phy_all$Concept_Name == "Systolic",
]
htn$Date <- as.Date(htn$Date)
htn_date_check <- htn[
  is.na(htn$htn_Date) | htn$htn_Date < htn$Date,
]

Dia <- htn_date_check[htn_date_check$Concept_Name == "Diastolic", ]
Sys <- htn_date_check[htn_date_check$Concept_Name != "Diastolic", ]

include<- "(Diastolic|Systolic)"

# Filter
phy_all <- phy_all[
  grepl(include, phy_all$Concept_Name, ignore.case = TRUE),]

wide_phy_all <- phy_all %>%
  pivot_wider(
    id_cols = EMPI,
    names_from = Concept_Name,
    values_from = Result
  )

fwrite(wide_phy_all,"wide_phy_all.csv")


wide_phy_all$Dia <- ifelse(
  wide_phy_all$EMPI %in% Dia$EMPI,
  0, 
  ifelse(wide_phy_all$Diastolic<90,0,1))

wide_phy_all$Sys <- ifelse(
  wide_phy_all$EMPI %in% Sys$EMPI,
  0, 
  ifelse(wide_phy_all$Systolic<140,0,1))

wide_phy_all<-merge(wide_phy_all,smok,by="EMPI",all.x=T)

wide_phy_all$keep<-ifelse(wide_phy_all$Dia==1 | wide_phy_all$Sys==1 | wide_phy_all$Smok==1,0,1 )
phy_remove<-wide_phy_all[wide_phy_all$keep==0,]
IDs$phy<-ifelse(IDs$EMPI %in% phy_remove$EMPI,1,0)
fwrite(wide_phy_all,"wide_phy_all.csv")



#################################################

rm(list = setdiff(ls(), "IDs"))

CAD_ID_beforemed$onset_age<-as.numeric(difftime(CAD_ID_beforemed$cad_Date, CAD_ID_beforemed$Date_of_Birth, units = "days")) / 365.25

IDs_med<-IDs[(IDs$phy==0 & IDs$lab==0),]
med_1<-fread("PN16_20250218_110101-6_Med.txt")
med_1<- med_1[!is.na(match(med_1$EMPI, IDs_med$EMPI))]

include_htn <- "(Amlodipine|Norvasc|Amiloride|Midamor|Atenolol|Tenormin|Acebutolol|Sectral|Azilsartan|Edarbi|
Aliskiren|Tekturna|Benazepril|Lotensin|Betaxolol|Kerlone|Bisoprolol|Zebeta|Bendroflumethiazide|
Bumetanide|Bumex|Candesartan|Atacand|Carvedilol|Coreg|Captopril|Capoten|Chlorthalidone|Hygroton|
Chlorothiazide|Diuril|Clonidine|Catapres|Diltiazem|Cardizem|Cartia|Tiazac|Doxazosin|Cardura|
Enalapril|Vasotec|Eplerenone|Inspra|Eprosartan|Teveten|Ethacrynic|Edecrin|Fosinopril|Monopril|
Felodipine|Plendil|Furosemide|Lasix|Guanfacine|Tenex|Intuniv|Hydralazine|Apresoline|
Hydrochlorothiazide|HCTZ|Microzide|Indapamide|Lozol|Isradipine|DynaCirc|Irbesartan|Avapro|
Labetalol|Normodyne|Trandate|Lisinopril|Prinivil|Zestril|Losartan|Cozaar|Metoprolol|Lopressor|Toprol|
Methyldopa|Aldomet|Moexipril|Univasc|Minoxidil|Loniten|Nadolol|Corgard|Nebivolol|Bystolic|
Nicardipine|Cardene|Nifedipine|Procardia|Adalat|Nisoldipine|Sular|Olmesartan|Benicar|
Pindolol|Visken|Prazosin|Minipress|Perindopril|Aceon|Propranolol|Inderal|
Quinapril|Accupril|Ramipril|Altace|Spironolactone|Aldactone|Telmisartan|Micardis|
Triamterene|Dyrenium|Timolol|Blocadren|Trandolapril|Mavik|Terazosin|Hytrin|
Torsemide|Demadex|Valsartan|Diovan|Verapamil|Calan|Isoptin|Verelan)"

include_lipid<-"(Alirocumab|Praluent|
Atorvastatin|Lipitor|Bempedoic|Nexletol|Cholestyramine|Questran|Colesevelam|Welchol|
Colestipol|Colestid|Evolocumab|Repatha|Ezetimibe|Zetia|Fenofibrate|Tricor|Triglide|Fluvastatin|Lescol|
Gemfibrozil|Lopid|Lovastatin|Mevacor|Pitavastatin|Livalo|Pravastatin|Pravachol|Niacin|Niaspan|
Lovaza|Omacor|Vascepa|Rosuvastatin|Crestor|Simvastatin|Zocor)"

med_1 <- med_1[
  grepl(paste(include_htn, include_lipid, sep = "|"), Medication, ignore.case = TRUE)
]

med_1$Medication[grepl(include_htn, med_1$Medication, ignore.case = TRUE)] <- "htn"
med_1$Medication[grepl(include_lipid, med_1$Medication, ignore.case = TRUE)] <- "lipid"


cad_dates<-IDs_med[,c(1,4)]
med_1 <- merge(med_1, cad_dates, by = "EMPI", all.x = TRUE)
setDT(med_1)
med_1[, Medication_Date := as.Date(Medication_Date, format = "%m/%d/%Y")]
med_1[, cad_Date := as.Date(cad_Date)]

# Compute distance
med_1[, diff_days := as.numeric(difftime(Medication_Date, cad_Date, units = "days"))]

med_1 <- med_1[
  ,
  {
    before <- .SD[diff_days < 0]
    
    if (nrow(before) > 0) {
      # choose closest BEFORE (largest negative)
      before[which.max(diff_days)]
    } else {
      after <- .SD[diff_days > 0]
      if (nrow(after) > 0) {
        # choose closest AFTER (smallest positive)
        after[which.min(diff_days)]
      } else {
        NULL
      }
    }
  },
  by = .(EMPI,Medication)
]

med_1<-med_1[,c(1,2,6,17,18)]

med_all<-rbind(med_1,med_2,med_3,med_4,med_5,med_6)

htn <- med_all[med_all$Medication == "htn",]
lipid<-med_all[med_all$Medication == "lipid",]

htn_date<-MGBB_Phenos_2CODE_2022_06_21[!is.na(match(MGBB_Phenos_2CODE_2022_06_21$EMPI, htn$EMPI)),]
lipid_date<-MGBB_Phenos_2CODE_2022_06_21[!is.na(match(MGBB_Phenos_2CODE_2022_06_21$EMPI, lipid$EMPI)),]
htn_date<-htn_date[,c(2,12)]
lipid_date<-lipid_date[,c(2,37)]
htn<-merge(htn,htn_date,by="EMPI")
lipid<-merge(lipid,lipid_date,by="EMPI")

htn_keep <- htn[
  is.na(htn$htn_Date) | htn$htn_Date > htn$Medication_Date,
]
lipid_keep <- lipid[
  is.na(lipid$hyperlipidemia_Date) | lipid$hyperlipidemia_Date > lipid$Medication_Date,
]
htn_keep$htn_Date<-NULL
lipid_keep$hyperlipidemia_Date<-NULL
med_all_final<-rbind(htn_keep,lipid_keep)

med<-as.data.frame(table(med_1$Medication))

CAD_ID_beforemed$onset_age<-as.numeric(difftime(CAD_ID_beforemed$cad_Date, CAD_ID_beforemed$Date_of_Birth, units = "days")) / 365.25

smurfless <- IDs_med[is.na(match(IDs_med$EMPI, med_all_final$EMPI)), ]
smurfless$onset<-as.numeric(difftime(smurfless$cad_Date, smurfless$Date_of_Birth, units = "days")) / 365.25
smurfless_early<-smurfless[smurfless$onset<55,]
fwrite(smurfless,"smurfless.csv")
fwrite(IDs,"IDs.csv")
race<-mgbb_smurfless_pheno_20250707[,c("EMPI","race_ethnicity")]
smurfless_early<-merge(smurfless_early,race,by="EMPI",all.x=T)
fwrite(smurfless_early,"smurfless_early.csv")
table(smurfless_early$race_ethnicity)
