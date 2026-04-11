
library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
HR2_Final <- read_excel("HR2_Final.xlsx", sheet = 1)  # or specify sheet name

names(HR2_Final)[names(HR2_Final) == "Exposure"] <- "PRS"

#HR Plot
# Get the PRS order based on Total group
prs_order <- HR2_Final %>%
  filter(Group == "All Ancestries") %>%
  arrange(desc(HR)) %>%
  pull(PRS)

# Apply the same PRS order to all groups
HR2_Final <- HR2_Final %>%
  mutate(PRS_ordered = factor(PRS, levels = prs_order))

# Plot
png("HR_plot.png", width = 2000, height = 1200, res = 150)
ggplot(HR2_Final, aes(x = PRS_ordered, y = HR, color = Group, shape = Group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c("British White Ancestry" = "steelblue",
                                "Non-British White Ancestry" = "orange",
                                "All Ancestries" = "red")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_rect(fill = "gray90"),
        legend.title = element_blank()) +
  labs(y = "HR", x = "")
dev.off()


#HR2 slop plot
library(ggplot2)
library(tidyr)
library(dplyr)
PRS_summary <- read_excel("PRS_summary.xlsx", sheet = 1)  # or specify sheet name

ranking<-prsnewsum[,c(1,3,5)]
colnames(ranking) <- c("PRS","Population","Individual")

# Calculate slope direction
ranking <- ranking %>%
  mutate(slope_dir = case_when(
    Individual < Population ~ "Positive",
    Individual > Population ~ "Negative",
    TRUE ~ "Constant"
  ))

# Convert data to long format for ggplot
ranking_long <- ranking %>%
  pivot_longer(cols = c("Population", "Individual"),
               names_to = "Rank_Type", values_to = "Rank")

# Plot slope graph
png("myplot_slope_p.png", width =8, height =12, units = "in", res = 300)
ggplot(ranking_long, aes(x = Rank_Type, y = Rank, group = PRS)) +
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
  geom_text(data = ranking_long %>% filter(Rank_Type == "Individual"),
            aes(label = PRS, x = 1), hjust = 1.1, size = 3) +
  geom_text(data = ranking_long %>% filter(Rank_Type == "Population"),
            aes(label = PRS, x = 2), hjust = -0.1, size =3) +
  coord_cartesian(clip = 'off')  
dev.off()

#HR2 slope plot
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
PRS_summary <- read_excel("PRS_summary.xlsx", sheet = 1)  # or specify sheet name

ranking_HR2<-PRS_summary[,c(1,2,3,4,5)]
colnames(ranking_HR2) <- c("PRS","HR","population","percentage","individual")

# Calculate slope direction
ranking_HR2 <- ranking_HR2 %>%
  mutate(slope_dir = case_when(
    individual < population ~ "Positive",
    individual > population ~ "Negative",
    TRUE ~ "Constant"
  ))

# Convert data to long format for ggplot
ranking_long_HR2 <- ranking_HR2 %>%
  pivot_longer(cols = c("population", "individual"),
               names_to = "Rank_Type", values_to = "Rank")

# Plot slope graph
png("slope_better.png", width = 1800, height = 2000, res = 200)
ggplot(ranking_long_HR2, aes(x = Rank_Type, y = Rank, group = PRS)) +
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
  geom_text(data = ranking_long_HR2 %>% filter(Rank_Type == "individual"),
            aes(label = PRS, x = 1), hjust = 1.1, size = 3) +
  geom_text(data = ranking_long_HR2 %>% filter(Rank_Type == "population"),
            aes(label = PRS, x = 2), hjust = -0.1, size =3) +
  coord_cartesian(clip = 'off')  
dev.off()
