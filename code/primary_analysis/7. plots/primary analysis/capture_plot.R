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
  filename = "my_plot.png",
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
install.packages("ggh4x")
library(ggh4x)

ggplot(summary_df_55_final,
       aes(x = PRS, y = prop, color = Group)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(
    aes(ymin = low_ci, ymax = high_ci),
    width = 0.6,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~ Group, scales = "free_x", ncol = 1) +
  facetted_pos_scales(
    x = list(
      NULL,   # top facet
      NULL,   # middle facet
      scale_x_discrete()  # bottom facet
    )
  ) +
  scale_color_manual(values = c(
    "Female" = "steelblue",
    "Male"   = "orange",
    "Total"  = "red"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(y = "Percentage (%)", x = "")



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

###########################################with CI
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)


summary_df_55_final<-results
colnames(summary_df_55_final) <-
  c("PRS", "Group", "OR", "low_ci", "high_ci")

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


total_order <- summary_df_55_final %>%
  filter(Group == "All Ancestries") %>%
  arrange(desc(OR)) %>%
  pull(PRS)

summary_df_55_final <- summary_df_55_final %>%
  mutate(PRS = factor(PRS, levels = total_order))


ggplot(summary_df_55_final,
       aes(x = PRS, y = OR, color = Group)) +
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
    y = "OR/SD",
    x = ""
  )


ggsave(
  filename = "my_plot_OR.png",
  plot = last_plot(),
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)

####################################################
library(ggplot2)
library(dplyr)

library(ggplot2)

plot_one_group <- function(data, group_name, color,
                           show_x = FALSE, show_y_label = FALSE) {
  
  ggplot(data, aes(x = PRS, y = prop)) +
    geom_point(size = 3, color = color) +
    geom_errorbar(
      aes(ymin = low_ci, ymax = high_ci),
      width = 0.6,
      color = color
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      
      # X-axis control
      axis.text.x  = if (show_x) element_text(angle = 90, hjust = 1) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.title.x = element_blank(),
      
      # Y-axis: always show ticks & text
      axis.text.y  = element_text(),
      axis.ticks.y = element_line(),
      axis.title.y = if (show_y_label) element_text() else element_blank(),
      
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = group_name,
      y = if (show_y_label) "Percentage (%)" else NULL
    )
}


library(dplyr)
library(patchwork)

p_total <- summary_df_55_final %>%
  filter(Group == "Total") %>%
  plot_one_group(
    group_name = "Total",
    color = "red",
    show_x = FALSE,
    show_y = FALSE
  )

p_female <- summary_df_55_final %>%
  filter(Group == "Female") %>%
  plot_one_group(
    group_name = "Female",
    color = "steelblue",
    show_x = FALSE,
    show_y = TRUE
  )

p_male <- summary_df_55_final %>%
  filter(Group == "Male") %>%
  plot_one_group(
    group_name = "Male",
    color = "orange",
    show_x = TRUE,
    show_y = FALSE
  )
final_plot <- p_total / p_female / p_male +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

final_plot
ggsave(
  filename = "plot_sex_new.png",
  plot = last_plot(),
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)
###################################################
plot_one_group <- function(data, group_name, color,
                           show_x = FALSE, show_y_label = FALSE) {
  
  ggplot(data, aes(x = PRS, y = prop)) +
    geom_point(size = 3, color = color) +
    geom_errorbar(
      aes(ymin = low_ci, ymax = high_ci),
      width = 0.6,
      color = color
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      
      # X-axis control
      axis.text.x  = if (show_x) element_text(angle = 90, hjust = 1) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.title.x = element_blank(),
      
      # Y-axis: always show ticks & text
      axis.text.y  = element_text(),
      axis.ticks.y = element_line(),
      axis.title.y = if (show_y_label) element_text() else element_blank(),
      
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = group_name,
      y = if (show_y_label) "Percentage (%)" else NULL
    )
}


library(dplyr)
library(patchwork)

p_total <- summary_df_55_final %>%
  filter(Group == "All Ancestries") %>%
  plot_one_group(
    group_name = "All Ancestries",
    color = "red",
    show_x = FALSE,
    show_y = FALSE
  )

p_white <- summary_df_55_final %>%
  filter(Group == "British White Ancestry") %>%
  plot_one_group(
    group_name = "British White Ancestry",
    color = "steelblue",
    show_x = FALSE,
    show_y = TRUE
  )

p_other <- summary_df_55_final %>%
  filter(Group == "Non-British White Ancestry") %>%
  plot_one_group(
    group_name = "Non-British White Ancestry",
    color = "orange",
    show_x = TRUE,
    show_y = FALSE
  )
final_plot <- p_total / p_white / p_other +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

final_plot
ggsave(
  filename = "plot_new.png",
  plot = last_plot(),
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)
