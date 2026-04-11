###########################################with CI
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)


summary_df_55_final<-boot_summary_df_55
colnames(summary_df_55_final) <-
  c("PRS", "Group", "prop", "low_ci", "high_ci")

summary_df_55_final <- summary_df_55_final %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Total", "British_White", "other")
    )
  ) %>%
  mutate(
    Group = case_when(
      Group == "Total" ~ "All Ancestries",
      Group == "British_White" ~ "British White Ancestry",
      Group == "other" ~ "Non-British White Ancestry",
      TRUE ~ as.character(Group)
    )
  )


summary_df_55_final <- summary_df_55_final %>%
  mutate(across(c(prop, low_ci, high_ci), ~ .x * 100))
total_order <- summary_df_55_final %>%
  filter(Group == "All Ancestries") %>%
  arrange(desc(prop)) %>%
  pull(PRS)

summary_df_55_final <- summary_df_55_final %>%
  mutate(PRS = factor(PRS, levels = total_order))


ggplot(summary_df_55_final,
       aes(x = PRS, y = prop, color = Group)) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = low_ci, ymax = high_ci),
    width = 0.6,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c(
    "British White Ancestry" = "steelblue",
    "Non-British White Ancestry" = "orange",
    "All Ancestries" = "red"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    
    # ✅ white backgrounds (fixes ggsave grey issue)
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    
    # ✅ KEEP grey facet strip
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text       = element_text(face = "bold"),
    
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = "Percentage (%)",
    x = ""
  )


ggsave(
  filename = "my_plot_p.png",
  plot = last_plot(),
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)
###########################################################


library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)


summary_df_55_final<-boot_summary_df_55_sex
colnames(summary_df_55_final) <-
  c("PRS", "Group", "prop", "low_ci", "high_ci")

summary_df_55_final <- summary_df_55_final %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Total", "Male", "Female")
    )
  ) 


summary_df_55_final <- summary_df_55_final %>%
  mutate(across(c(prop, low_ci, high_ci), ~ .x * 100))
total_order <- summary_df_55_final %>%
  filter(Group == "Total") %>%
  arrange(desc(prop)) %>%
  pull(PRS)

summary_df_55_final <- summary_df_55_final %>%
  mutate(PRS = factor(PRS, levels = total_order))

ggplot(summary_df_55_final,
       aes(x = PRS, y = prop, color = Group)) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = low_ci, ymax = high_ci),
    width = 0.6,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c(
    "Female" = "steelblue",
    "Male"   = "orange",
    "Total"  = "red"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = "Percentage (%)",
    x = ""
  )


ggsave(
  filename = "my_plot_sex.png",
  plot = last_plot(),
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)



