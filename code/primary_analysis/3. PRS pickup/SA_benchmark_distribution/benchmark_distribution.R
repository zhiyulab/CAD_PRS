# Quantile edges for top 10%, 10–20%, ..., 90–100%
probs <- seq(0, 1, by = 0.1)  # 0.0, 0.1, 0.2, ..., 1.0

common_counts <- numeric(length(probs)-1)
bin_labels <- character(length(probs)-1)

for (i in 1:(length(probs)-1)) {
  lower <- quantile(UKB_PRS$df53, probs[i], na.rm = TRUE)
  upper <- quantile(UKB_PRS$df53, probs[i+1], na.rm = TRUE)
  
  # IDs in this bin
  bin_ids <- UKB_PRS$id[UKB_PRS$df53 >= lower & UKB_PRS$df53 < upper]
  
  # count how many are also in 'add'
  common_counts[i] <- sum(bin_ids %in% cohort_p_55$id)
  
  # label for plotting
  bin_labels[i] <- paste0(probs[i]*100, "-", probs[i+1]*100, "%")
}

# create dataframe
df_plot <- data.frame(
  bin = factor(bin_labels, levels = rev(bin_labels)),  # reverse so top bin is first
  common_count = common_counts
)

library(ggplot2)

ggplot(df_plot, aes(x = bin, y = common_count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "PRS Percentile Bin", y = "Number of Common IDs",
       title = "Number of IDs in 'add' by PRS Percentile Bin") +
  theme_minimal() +
  coord_flip()  # optional: flip to horizontal for better readability