
library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggplot2)
library(tidyr)
library(dplyr)
prsnewsum <- read_excel("prsnewsum", sheet = 1)  # or specify sheet name

ranking<-prsnewsum[,c(2,4,13)]
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
library(ggplot2)
library(dplyr)

library(ggplot2)
library(dplyr)

p <- ggplot(ranking_long, aes(x = Rank_Type, y = Rank, group = PRS)) +
  geom_line(aes(color = slope_dir), size = 1.2) +
  geom_point(aes(color = slope_dir), size = 4) +
  scale_y_reverse() +  # reverse y-axis so rank 1 is top
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
    legend.position = "top",
    plot.background = element_rect(fill = "white", color = NA)  # white background
  ) +
  geom_text(
    data = ranking_long %>% filter(Rank_Type == "Individual"),
    aes(label = PRS, x = 1),
    hjust = 1.1,
    size = 3
  ) +
  geom_text(
    data = ranking_long %>% filter(Rank_Type == "Population"),
    aes(label = PRS, x = 2),
    hjust = -0.1,
    size = 3
  ) +
  coord_cartesian(clip = 'off')

# Save the plot with white background
ggsave("myplot_slope_p.png", plot = p, width = 8, height = 12, units = "in", dpi = 300, bg = "white")
