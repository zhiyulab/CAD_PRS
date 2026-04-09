#########################################################
UKB_PRS<-UKB_PRS[,c("id","df53")]
cutoff <- quantile(UKB_PRS$df53, 0.9, na.rm = TRUE)
UKB_PRS$cutoff_10 <- ifelse(UKB_PRS$df53 >= cutoff, 1, 0)

cutoff <- quantile(UKB_PRS$df53, 0.8, na.rm = TRUE)
UKB_PRS$cutoff_20 <- ifelse(UKB_PRS$df53 >= cutoff, 1, 0)

cutoff <- quantile(UKB_PRS$df53, 0.99, na.rm = TRUE)
UKB_PRS$cutoff_1 <- ifelse(UKB_PRS$df53 >= cutoff, 1, 0)

cutoff <- quantile(UKB_PRS$df53, 0.95, na.rm = TRUE)
UKB_PRS$cutoff_5 <- ifelse(UKB_PRS$df53 >= cutoff, 1, 0)

#precision-recall

library(dplyr)
library(tidyr)
library(ggplot2)

colnames(UKB_PRS)<-c("id","df53","10% cutoff","20% cutoff","1% cutoff","5% cutoff","case")
UKB_PRS$case<-ifelse(UKB_PRS$id %in% cohort_p_55$id,1,0)

threshold_cols <- c("1% cutoff", "5% cutoff", "10% cutoff", "20% cutoff")


pr_df <- UKB_PRS %>%
  pivot_longer(cols = all_of(threshold_cols),
               names_to = "threshold",
               values_to = "pred") %>%
  group_by(threshold) %>%
  summarise(
    TP = sum(pred == 1 & case == 1),
    FP = sum(pred == 1 & case == 0),
    FN = sum(pred == 0 & case == 1),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    .groups = "drop"
  )

library(ggplot2)
library(ggrepel)

pr_plot <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_line(color = "grey40") +
  geom_point(size = 3, color = "steelblue") +
  geom_text_repel(
    aes(label = threshold),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # white panel
    plot.background = element_rect(fill = "white", color = NA)   # white overall plot
  ) +
  labs(
    x = "Recall",
    y = "Precision"
  )

ggsave(
  filename = "precision_recall_curve.png",
  plot = pr_plot,
  width = 6,
  height = 4,
  dpi = 300
)