library(data.table)
pheno<-fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt") 
pheno = data.frame(pheno)
pheno$IDs_toRemove_SampleQCv3 = ifelse((!(pheno$Submitted_Gender == pheno$Inferred_Gender) | pheno$Non_Consented== 1),1,0) 
pheno = pheno[which(pheno$IDs_toRemove_SampleQCv3==0),]

#get BMI column
bmi<-pheno[,c(1,18)]
library(readr)
SMuRFless_55 <- read_csv("/medpop/esp2/sliang/SMuRFless_CAD_under55.csv")
SMuRFless_55<-SMuRFless_55[,1]

#take top 10% cutoff as before with PRS
cutoff <- quantile(bmi$BMI, 0.9, na.rm = TRUE)
bmi$cutoff_10 <- ifelse(bmi$BMI >= cutoff, 1, 0)
SMuRFless_55_BMI<-merge(SMuRFless_55,bmi,by="id",all.x=TRUE)
print(sum(SMuRFless_55_BMI$cutoff_10==1,na.rm=T))
print(sum(is.na(SMuRFless_55_BMI$cutoff_10)))

#get lpa and crp colulmns
lpa_crp<-pheno[,c(1,180,188)]
colnames(lpa_crp)<-c("id","crp","lpa")

#take top 10% cutoff as before with PRS
lpa_cutoff <- quantile(lpa_crp$lpa, 0.9, na.rm = TRUE)
crp_cutoff<-quantile(lpa_crp$crp, 0.9, na.rm = TRUE)
lpa_crp$lpa_cutoff<-ifelse(lpa_crp$lpa >= lpa_cutoff, 1, 0)
lpa_crp$crp_cutoff<-ifelse(lpa_crp$crp >= crp_cutoff, 1, 0)
SMuRFless_55_lpa_crp<-merge(SMuRFless_55,lpa_crp,by="id",all.x=TRUE)
print(sum(SMuRFless_55_lpa_crp$lpa_cutoff==1,na.rm=T))
print(sum(is.na(SMuRFless_55_lpa_crp$lpa_cutoff)))
print(sum(SMuRFless_55_lpa_crp$crp_cutoff==1,na.rm=T))
print(sum(is.na(SMuRFless_55_lpa_crp$crp_cutoff)))

#make a stacked bar plot
library(ggplot2)
library(tidyr)

#percentage already calcaulted as (number of cutoff==1) / (size of SMuRFless early onset group - number of no score)
df <- data.frame(
  RiskFactor = c("BMI", "Lp(a)", "CRP", "PRS"),
  Capture = c(12.1, 12.9, 12.2, 31.8)
)

# Add Remaining
df$Remaining <- 100 - df$Capture

# Convert to long format
df_long <- pivot_longer(df, cols = c("Capture", "Remaining"), names_to = "Type", values_to = "Value")

# Reorder to make Capture come first (so it's on the left)
df_long$Type <- factor(df_long$Type, levels = c("Remaining", "Capture"))

png("myplot.png", width = 6, height = 4, units = "in", res = 300)

ggplot(df_long, aes(x = RiskFactor, y = Value, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("Remaining" = "grey90", "Capture" = "steelblue")) +
  labs(
    x = "",
    y = "Percentage (%) of early onset unexpected CAD captured in top decile"
  ) +
  geom_text(
    data = subset(df_long, Type == "Capture"),
    aes(label = paste0(Value, "%"), y = Value / 2),
    color = "white",
    size = 3
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),  # remove minor grid lines
    panel.border = element_blank(),      # remove outer border (optional)
    legend.position = "none"
  )


dev.off()  

