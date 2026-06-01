# Script for Accumulation Analysis of FHV RNA Results

# Virus: Flock House Virus (FHV) [RefSeq: NC_004144.1]
# Titre: 10E8 IU/mL
# Fly: w1118 Control vs miRNA Mutant
# Collection Time: 3-dpi
# Amplicon: 139 bp, targeting [RefSeq: NP_689442.1]
# Fig. 5B, Table S5

# Rationale
#  Determine whether miRNA mutant line changes FHV viral RNA accumulation compared to paired control.

# Approach
#  Differences in FHV viral RNA accumulation between miRNA mutant and control were analysed using Mean
#  Normalised Expression (MNE) and Paired T-Tests. MNE was calculated using the Q-Gene formula (1).
#  For each pairwise comparison between miRNA mutant and control, a Paired Two-tailed T-test was used.
#  Due to the targeted nature of this experiment, we did not adjust for multiple comparisons, unlike
#  the screen (2). Confidence intervals for Log2FC were estimated by bootstrapping.

# References
# 1. Simon P. 2003. Q-Gene: processing quantitative real-time RT–PCR data. Bioinformatics 19:1439–1440.
# 2. Hooper R. 2025. To adjust, or not to adjust, for multiple comparisons. J Clin Epidemiol. 180:111688.

# Preparation
library(dplyr)
library(purrr)
library(ggplot2)

# Set seed for reproducibility of bootstrap sampling
set.seed(123)

# Load MNE data from CSV and ensure WT/MT are numeric
df <- read.csv("FHV/DATA.csv") %>%
  mutate(
    WT = as.numeric(WT),
    MT = as.numeric(MT)
  )

# Number of bootstrap resamples
B = 10000

# Function to compute bootstrap distribution of log2 fold change
bootstrap_log2fc <- function(wt, mt, B = B) {
  
  # Number of paired observations
  n <- length(wt)
  
  # Vector to store bootstrap estimates
  boot_vals <- numeric(B)
  
  # Bootstrap resampling loop
  for (i in seq_len(B)) {
    
    # Sample indices with replacement (paired bootstrap)
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    
    # Resampled WT and MT values
    wt_b <- wt[idx]
    mt_b <- mt[idx]
    
    # Compute log2 fold change for bootstrap sample
    boot_vals[i] <- log2(mean(mt_b, na.rm = TRUE) /
                           mean(wt_b, na.rm = TRUE))
  }
  
  # Return full bootstrap distribution
  return(boot_vals)
}

# Main analysis: compute statistics per miRNA
final_results <- df %>%
  
  # Group data by miRNA
  group_by(miRNA) %>%
  
  summarise(
    
    # Mean expression in each group
    WT_mean = mean(WT, na.rm = TRUE),
    MT_mean = mean(MT, na.rm = TRUE),
    
    # Point estimate of log2 fold change
    log2FC = log2(MT_mean / WT_mean),
    
    # Paired t-test between WT and MT
    test = list(t.test(WT, MT, paired = TRUE)),
    
    # Bootstrap distribution of log2FC
    boot = list(bootstrap_log2fc(WT, MT, B = B)),
    
    .groups = "drop"
  ) %>%
  
  # Extract stats
  mutate(
    p_value = sapply(test, function(x) x$p.value),
    statistic = sapply(test, function(x) unname(x$statistic)),
    
    # Bootstrap 95% confidence intervals
    ci_low = sapply(boot, function(x) quantile(x, 0.025, na.rm = TRUE)),
    ci_high = sapply(boot, function(x) quantile(x, 0.975, na.rm = TRUE))
  ) %>%
  
  # Keep only final output columns
  select(
    miRNA,
    WT_mean,
    MT_mean,
    log2FC,
    ci_low,
    ci_high,
    statistic,
    p_value
  )


# Save results to CSV file
write.csv(
  final_results,
  "Table S5. LOG2FC_FHV.csv",
  row.names = FALSE
)

# Print final results
print(final_results)

# Format miRNA as factor
final_results$miRNA <- as.factor(final_results$miRNA)

# Create significance column from adjusted p-value
final_results$signif <- ifelse(
  final_results$p_value <= 0.05,
  "Significant",
  "Not Significant"
)

# Order miRNAs by log2FC
ordered_miRNAs <- as.character(
  final_results$miRNA[order(final_results$log2FC)]
)

# Set plotting order
final_results$miRNA <- factor(
  as.character(final_results$miRNA),
  levels = unique(ordered_miRNAs)
)

# Reformat significance labels for plot
final_results$signif <- ifelse(
  final_results$signif == "Significant",
  "Significant (p ≤ 0.05)",
  "Not Significant (p > 0.05)"
)

# Plot limits
x_limits <- c(-8, 8)
x_breaks <- seq(-8, 8, by = 1)

# Plot Fig. 5B
ggplot(final_results, aes(x = log2FC, y = miRNA, color = signif)) +
  
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.25, linewidth = 1.2
  ) +
  
  geom_point(size = 3.8, stroke = 1) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "gray40",
    linewidth = 1.1
  ) +
  
  scale_color_manual(values = c(
    "Significant (p ≤ 0.05)" = "#1f78b4",
    "Not Significant (p > 0.05)" = "#e31a1c"
  )) +
  
  scale_x_continuous(limits = x_limits, breaks = x_breaks) +
  scale_y_discrete(drop = FALSE) +
  
  labs(
    x = "FHV RNA Accumulation (Log2 Fold Change)",
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

# Save Fig. 5B
ggsave(
  filename = "Fig. 5B.tiff",
  plot = last_plot(),
  device = "tiff",
  units = "px",
  width = 5000,
  height = 3000,
  dpi = 300,
  compression = "lzw"
)
