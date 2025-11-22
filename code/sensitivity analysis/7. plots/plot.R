library(tidyr)
library(dplyr)


library(readxl)
summary_df_55_final <- read_excel("summary_df_55_final.xlsx", sheet = 1)  # or specify sheet name

long_df <- summary_df_55_final %>%
  pivot_longer(
    cols = c(Total, British_White, Non_British_White),
    names_to = "Group",
    values_to = "Percent"
  ) %>%
  mutate(
    Group = case_when(
      Group == "Total" ~ "All Ancestries",
      Group == "British_White" ~ "British White Ancestry",
      Group == "Non_British_White" ~ "Non-British White Ancestry",
      TRUE ~ Group
    )
  )
long_df$Percent<-long_df$Percent*100


library(dplyr)
library(ggplot2)

# Assuming your long_df has columns: PRS, Percent, Group

# Step 1: Get Total Percent per PRS
total_order <- long_df %>%
  filter(Group == "All Ancestries") %>%
  arrange(desc(Percent)) %>%
  select(PRS) %>%
  pull()

# Step 2: Make PRS an ordered factor in long_df based on total_order
long_df <- long_df %>%
  mutate(PRS = factor(PRS, levels = total_order))

# Step 3: Plot with PRS ordered by Total Percent
png("individual_plot.png", width = 2000, height = 1200, res = 150)
ggplot(long_df, aes(x = PRS, y = Percent, color = Group, shape = Group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c("British White Ancestry" = "steelblue",
                                "Non-British White Ancestry" = "orange",
                                "All Ancestries" = "red")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_rect(fill = "gray90"),
        legend.title = element_blank()) +
  labs(y = "Percentage", x = "")
# plot code here
dev.off()

total_order <- HR2_Final %>%
  filter(Group == "All Ancestries") %>%
  arrange(desc(HR)) %>%
  select(Exposure) %>%
  pull()

HR2_Final <- HR2_Final %>%
  mutate(Exposure = factor(Exposure, levels = total_order))

# Step 3: Plot with PRS ordered by Total Percent
png("myplot_HR.png", width =12, height = 8, units = "in", res = 300)
ggplot(HR2_Final, aes(x = Exposure, y = HR, color = Group, shape = Group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c("British White Ancestry" = "steelblue",
                                "Non-British White Ancestry" = "orange",
                                "All Ancestries" = "red")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_rect(fill = "gray90"),
        legend.title = element_blank()) +
  labs(y = "HR/SD", x = "")
# plot code here
dev.off()