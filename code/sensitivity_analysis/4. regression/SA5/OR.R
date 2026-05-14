library(data.table)
library(readr)
pheno<-fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt") 

regression<-mgbb_smurfless_pheno_20250707
genotype<-MGB_PRS_58[,(1:59)]
regression<-merge(regression,genotype,by = "EMPI", all.y = TRUE)
regression<-regression[,c(1:4)]
rm(list = setdiff(ls(), c("regression","genotype")))


regression$outcome<-ifelse(regression$EMPI %in% IDs$EMPI,1,0)
regression$Date_of_Birth <- as.Date(regression$Date_of_Birth, format = "%m/%d/%Y")
regression$age <- as.Date("2026-01-17") - regression$Date_of_Birth
regression$age <- as.numeric(regression$age) / 365.25

regression_other <- regression[
  (is.na(regression$race_ethnicity) | (!is.na(regression$race_ethnicity) & regression$race_ethnicity != "Non-Hispanic White")),
]
regression_white<-regression[!is.na(regression$race_ethnicity) & (regression$race_ethnicity=="Non-Hispanic White"),]

##############################
library(tibble)        # For tidy data frames
library(dplyr)         # For data manipulation (if needed)
library(purrr)         # For functional programming tools (e.g., map)
library(furrr)         # For parallel versions of purrr functions
library(future)        # Backend for parallel computing

# Define your exposure variables (replace with actual names if not df1 to df58)
exposures <- paste0("df", 1:58)

# Covariates string (update as needed)
covariates <- "age + Gender_Legal_Sex + race_ethnicity"

run_logistic <- function(var) {
  formula <- as.formula(
    paste0("outcome ~ ", var, " + ", covariates)
  )
  
  model <- tryCatch(
    glm(formula, data = regression, family = binomial),
    error = function(e) return(NULL)
  )
  
  if (is.null(model)) {
    return(tibble(Exposure = var, Term = NA, OR = NA, CI_lower = NA, CI_upper = NA, P_value = NA))
  }
  
  coefs <- summary(model)$coefficients
  
  # Find the row that matches the exposure variable name
  if (!(var %in% rownames(coefs))) {
    return(tibble(Exposure = var, Term = NA, OR = NA, CI_lower = NA, CI_upper = NA, P_value = NA))
  }
  
  row <- coefs[var, ]
  or <- exp(row["Estimate"])
  se <- row["Std. Error"]
  ci_low <- exp(row["Estimate"] - 1.96 * se)
  ci_high <- exp(row["Estimate"] + 1.96 * se)
  p <- row["Pr(>|z|)"]
  
  tibble(
    Exposure = var,
    Term = var,
    OR = or,
    CI_lower = ci_low,
    CI_upper = ci_high,
    P_value = p
  )
}

# Run all models in parallel
logistic_results <- future_map_dfr(exposures, run_logistic, .progress = TRUE)

##############################
library(tibble)        # For tidy data frames
library(dplyr)         # For data manipulation (if needed)
library(purrr)         # For functional programming tools (e.g., map)
library(furrr)         # For parallel versions of purrr functions
library(future)        # Backend for parallel computing

# Define your exposure variables (replace with actual names if not df1 to df58)
exposures <- paste0("df", 1:58)

# Covariates string (update as needed)
covariates <- "age + Gender_Legal_Sex + race_ethnicity"

run_logistic <- function(var) {
  formula <- as.formula(
    paste0("outcome ~ ", var, " + ", covariates)
  )
  
  model <- tryCatch(
    glm(formula, data = regression_other, family = binomial),
    error = function(e) return(NULL)
  )
  
  if (is.null(model)) {
    return(tibble(Exposure = var, Term = NA, OR = NA, CI_lower = NA, CI_upper = NA, P_value = NA))
  }
  
  coefs <- summary(model)$coefficients
  
  # Find the row that matches the exposure variable name
  if (!(var %in% rownames(coefs))) {
    return(tibble(Exposure = var, Term = NA, OR = NA, CI_lower = NA, CI_upper = NA, P_value = NA))
  }
  
  row <- coefs[var, ]
  or <- exp(row["Estimate"])
  se <- row["Std. Error"]
  ci_low <- exp(row["Estimate"] - 1.96 * se)
  ci_high <- exp(row["Estimate"] + 1.96 * se)
  p <- row["Pr(>|z|)"]
  
  tibble(
    Exposure = var,
    Term = var,
    OR = or,
    CI_lower = ci_low,
    CI_upper = ci_high,
    P_value = p
  )
}

# Run all models in parallel
logistic_results_other <- future_map_dfr(exposures, run_logistic, .progress = TRUE)

##############################
library(tibble)        # For tidy data frames
library(dplyr)         # For data manipulation (if needed)
library(purrr)         # For functional programming tools (e.g., map)
library(furrr)         # For parallel versions of purrr functions
library(future)        # Backend for parallel computing

# Define your exposure variables (replace with actual names if not df1 to df58)
exposures <- paste0("df", 1:58)

# Covariates string (update as needed)
covariates <- "age + Gender_Legal_Sex"

run_logistic <- function(var) {
  formula <- as.formula(
    paste0("outcome ~ ", var, " + ", covariates)
  )
  
  model <- tryCatch(
    glm(formula, data = regression_white, family = binomial),
    error = function(e) return(NULL)
  )
  
  if (is.null(model)) {
    return(tibble(Exposure = var, Term = NA, OR = NA, CI_lower = NA, CI_upper = NA, P_value = NA))
  }
  
  coefs <- summary(model)$coefficients
  
  # Find the row that matches the exposure variable name
  if (!(var %in% rownames(coefs))) {
    return(tibble(Exposure = var, Term = NA, OR = NA, CI_lower = NA, CI_upper = NA, P_value = NA))
  }
  
  row <- coefs[var, ]
  or <- exp(row["Estimate"])
  se <- row["Std. Error"]
  ci_low <- exp(row["Estimate"] - 1.96 * se)
  ci_high <- exp(row["Estimate"] + 1.96 * se)
  p <- row["Pr(>|z|)"]
  
  tibble(
    Exposure = var,
    Term = var,
    OR = or,
    CI_lower = ci_low,
    CI_upper = ci_high,
    P_value = p
  )
}

# Run all models in parallel
logistic_results_white <- future_map_dfr(exposures, run_logistic, .progress = TRUE)

write.csv(logistic_results_white, "logistic_results_white.csv", row.names = FALSE)
write.csv(logistic_results_other, "logistic_results_other.csv", row.names = FALSE)
write.csv(logistic_results, "logistic_results.csv", row.names = FALSE)

##############################
library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
colnames(OR_final)[1] <- "PRS"

#HR Plot
# Get the PRS order based on Total group
prs_order <- OR_final %>%
  filter(Group == "Total") %>%
  arrange(desc(OR)) %>%
  pull(PRS)

# Apply the same PRS order to all groups
OR_final <- OR_final %>%
  mutate(PRS_ordered = factor(PRS, levels = prs_order))

# Plot
png("HR_plot.png", width = 2000, height = 1200, res = 150)
ggplot(OR_final, aes(x = PRS_ordered, y = OR, color = Group, shape = Group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c("Subset_White" = "steelblue",
                                "Subset_Other" = "orange",
                                "Total" = "red")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_rect(fill = "gray90"),
        legend.title = element_blank()) +
  labs(y = "OR", x = "")
dev.off()


#HR2 slop plot
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
PRS_summary <- read_excel("PRS_summary.xlsx", sheet = 1)  # or specify sheet name

spearman_corr <- cor(PRS_summary$MGB_percentage_ranking_white, PRS_summary$UKB_percentage_ranking, method = "spearman")
print(spearman_corr)

ranking<-PRS_summary[,c(1,2,4)]
colnames(ranking) <- c("PRS","population","individual")

# Calculate slope direction
ranking <- ranking %>%
  mutate(slope_dir = case_when(
    individual < population ~ "Positive",
    individual > population ~ "Negative",
    TRUE ~ "Constant"
  ))

# Convert data to long format for ggplot
ranking <- ranking %>%
  pivot_longer(cols = c("individual", "population"),
               names_to = "Rank_Type", values_to = "Rank")

# Plot slope graph
png("slope_better.png", width = 1800, height = 2000, res = 200)
ggplot(ranking, aes(x = Rank_Type, y = Rank, group = PRS)) +
  geom_line(aes(color = slope_dir), size = 1.2) +
  geom_point(aes(color = slope_dir), size = 4) +
  scale_y_reverse() +  # reverse y-axis so rank 1 is top (keep for spacing)
  scale_color_manual(values = c(
    "Positive" = "red",
    "Negative" = "blue",
    "Constant" = "gray"
  )) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  geom_text(data = ranking %>% filter(Rank_Type == "individual"),
            aes(label = PRS, x = 1), hjust = 1.1, size = 3) +
  geom_text(data = ranking %>% filter(Rank_Type == "population"),
            aes(label = PRS, x = 2), hjust = -0.1, size =3) +
  coord_cartesian(clip = 'off')  
dev.off()


