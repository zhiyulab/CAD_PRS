library(data.table)
pheno<-fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt") 
pheno = data.frame(pheno)
pheno$IDs_toRemove_SampleQCv3 = ifelse((!(pheno$Submitted_Gender == pheno$Inferred_Gender) | pheno$Non_Consented== 1),1,0) 
pheno = pheno[which(pheno$IDs_toRemove_SampleQCv3==0),]

#get BMI column
bmi<-pheno[,c(1,18)]
colnames(bmi)
library(readr)
SMuRFless_55 <- read.csv("/home/unix/sliang/sliang/smurfless/SMuRFless_CAD_under55.csv")

#take top 10% cutoff as before with PRS
cutoff <- quantile(bmi$BMI, 0.9, na.rm = TRUE)
bmi$cutoff_10 <- ifelse(bmi$BMI >= cutoff, 1, 0)
SMuRFless_55_BMI<-merge(SMuRFless_55,bmi,by="id",all.x=TRUE)

###################################################
library(boot)

# Select the column
x <- SMuRFless_55_BMI$cutoff_10
x <- x[!is.na(x)]  # remove NAs

# Run bootstrap
boot_out <- boot(
  data = x,
  statistic = function(x, indices) {
    mean(x[indices] == 1)  # proportion of 1s
  },
  R = 1000  # number of bootstrap resamples
)

# Compute percentile confidence interval
ci <- boot.ci(boot_out, type = "perc")

# Collect results in a data frame
results_bmi <- data.frame(
  variable = "bmi_cutoff",
  prop = mean(x == 1),
  lower_95 = ci$percent[4],
  upper_95 = ci$percent[5]
)

fwrite(results_bmi,"results_bmi.csv")

###############################################################

#get lpa and crp colulmns
lpa_crp<-pheno[,c(1,180,188)]
colnames(lpa_crp)
colnames(lpa_crp)<-c("id","crp","lpa")

#take top 10% cutoff as before with PRS
lpa_cutoff <- quantile(lpa_crp$lpa, 0.9, na.rm = TRUE)
crp_cutoff<-quantile(lpa_crp$crp, 0.9, na.rm = TRUE)
lpa_crp$lpa_cutoff<-ifelse(lpa_crp$lpa >= lpa_cutoff, 1, 0)
lpa_crp$crp_cutoff<-ifelse(lpa_crp$crp >= crp_cutoff, 1, 0)
SMuRFless_55_lpa_crp<-merge(SMuRFless_55,lpa_crp,by="id",all.x=TRUE)



###################################################
library(boot)

# Select the column
x <- SMuRFless_55_lpa_crp$lpa_cutoff
x <- x[!is.na(x)]  # remove NAs

# Run bootstrap
boot_out <- boot(
  data = x,
  statistic = function(x, indices) {
    mean(x[indices] == 1)  # proportion of 1s
  },
  R = 1000  # number of bootstrap resamples
)

# Compute percentile confidence interval
ci <- boot.ci(boot_out, type = "perc")

# Collect results in a data frame
results_lpa <- data.frame(
  variable = "lpa_cutoff",
  prop = mean(x == 1),
  lower_95 = ci$percent[4],
  upper_95 = ci$percent[5]
)

fwrite(results_lpa,"results_lpa.csv")

###############################################################
library(boot)

# Select the column
x <- SMuRFless_55_lpa_crp$crp_cutoff
x <- x[!is.na(x)]  # remove NAs

# Run bootstrap
boot_out <- boot(
  data = x,
  statistic = function(x, indices) {
    mean(x[indices] == 1)  # proportion of 1s
  },
  R = 1000  # number of bootstrap resamples
)

# Compute percentile confidence interval
ci <- boot.ci(boot_out, type = "perc")

# Collect results in a data frame
results_crp <- data.frame(
  variable = "crp_cutoff",
  prop = mean(x == 1),
  lower_95 = ci$percent[4],
  upper_95 = ci$percent[5]
)

fwrite(results_crp,"results_crp.csv")

stop()
########################################################
library(ggplot2)
library(dplyr)


# Example: rename variables in your df
results_bmi <- results_BMI %>%
  mutate(
    variable = recode(variable,
                      "prs_lpa_cutoff" = "PRS + Lp(a)",
                      "prs_pce_cutoff" = "PRS + PCE"))


# Create the plot and assign it to a variable
p <- ggplot(results_bmi) +
  # Grey background bar (total = 1)
  geom_bar(aes(x = variable, y = 1), stat = "identity", fill = "grey90") +
  # Colored proportion bar
  geom_bar(aes(x = variable, y = prop), stat = "identity", fill = "steelblue") +
  # Add confidence intervals
  geom_errorbar(aes(x = variable, ymin = lower_95, ymax = upper_95), width = 0.2) +
  # Add proportion numbers on top of bars
  geom_text(aes(x = variable, y = prop, label = round(prop, 2)), 
            hjust = -1.5, size = 4) +
  coord_flip(clip = "off") +
  theme_minimal(base_size = 12) +
  labs(y = "Proportion (%)", x = "Variable") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.background = element_rect(fill = "white"),       # white panel
    plot.background = element_rect(fill = "white"),        # white plot
    panel.grid.major = element_blank(),                    # remove major grid
    panel.grid.minor = element_blank()                     # remove minor grid
  )

# Save the plot
ggsave("barplot_with_CI.png", plot = p, width = 6, height = 4, units = "in", dpi = 300, bg = "white")
