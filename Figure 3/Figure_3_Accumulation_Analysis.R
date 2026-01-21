# Experimental details
#  Virus: SINV
#  Titre: 10E5 IU/mL
#  Fly: w1118 Control vs miRNA Mutant
#  Fig: 3

# Rationale
#  To determine whether miRNA mutant line changes SINV RNA accumulation compared to paired control.

# Approach
#  Log base 2 transform data (helps with distribution given the data is virus accumulation)
#  Use Two-tailed Paired T-Test as previously established in (10.1016/j.virol.2016.12.009).

# Preparation
library(dplyr)
library(purrr)
library(ggplot2)

# Importing Data
data <- read.csv("SINV/SINV.csv", stringsAsFactors = FALSE)

# Analysis
# Log2-transform WT and MT
data <- data %>%
  mutate(
    WT_log2 = log2(WT + 1),
    MT_log2 = log2(MT + 1)
  )

# Perform paired t-test for miRNA-control
results <- data %>%
  group_by(miRNA) %>%
  summarise(
    mean_WT_log2 = mean(WT_log2, na.rm = TRUE),
    mean_MT_log2 = mean(MT_log2, na.rm = TRUE),
    t_test = list(t.test(WT_log2, MT_log2, paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ as.numeric(.x$statistic)),
    p_value = map_dbl(t_test, ~ .x$p.value),
    conf_low_log2 = map_dbl(t_test, ~ .x$conf.int[1]),
    conf_high_log2 = map_dbl(t_test, ~ .x$conf.int[2]),
    log2_fold_change = mean_WT_log2 - mean_MT_log2
  ) %>%
  select(
    miRNA,
    mean_WT_log2,
    mean_MT_log2,
    t_statistic,
    p_value,
    log2_fold_change,
    conf_low_log2,
    conf_high_log2
  ) %>%
  arrange(log2_fold_change)
print(results)

# Export Results
write.csv(results, "stat_results_SINV.csv", row.names = FALSE)

# Next: Plot T-Test results for Fig 3
final_results <- results
# Create significance column from p-value
final_results$signif <- ifelse(
  final_results$p_value <= 0.05,
  "Significant",
  "Not Significant"
)

# Order based on log2 fold change
ordered_miRNAs <- as.character(
  final_results$miRNA[order(final_results$log2_fold_change)]
)

# Add extra miRNAs at the bottom
# These are the two miRNAs which had no detectable SINV signal, and therefore no plottable data points. Adding them to bottom of Fig. 3.
extra_miRNAs <- c("miR-10", "miR-957")

final_results$miRNA <- factor(
  as.character(final_results$miRNA),
  levels = unique(c(extra_miRNAs, ordered_miRNAs))
)

final_results$signif <- ifelse(
  final_results$signif == "Significant",
  "Significant (p ≤ 0.05)",
  "Not Significant (p > 0.05)"
)

# Specify plot limits
x_limits <- c(-2, 2)
x_breaks <- seq(-2, 2, by = 1)

# Plot Fig. 3
ggplot(final_results, aes(x = log2_fold_change, y = miRNA, color = signif)) +
  
  geom_errorbarh(
    aes(xmin = conf_low_log2, xmax = conf_high_log2),
    height = 0.25, linewidth = 1.2
  ) +
  
  geom_point(size = 3.8, stroke = 1) +
  
  geom_vline(
    xintercept = 0, linetype = "dotted",
    color = "gray40", linewidth = 1.1
  ) +
  
  scale_color_manual(values = c(
    "Significant (p ≤ 0.05)" = "#1f78b4",
    "Not Significant (p > 0.05)" = "#e31a1c"
  )) +
  
  scale_x_continuous(limits = x_limits, breaks = x_breaks) +
  scale_y_discrete(drop = FALSE) +
  
  labs(
    x = "Log2 Fold Change",
    y = "miRNA mutant",
    color = "Significance"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 15),
    legend.text  = element_text(face = "bold", size = 14)
  )
