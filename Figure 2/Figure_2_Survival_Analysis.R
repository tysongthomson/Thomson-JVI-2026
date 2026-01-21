# Experimental details
#  Virus: FHV
#  Titre: 10E8 IU/mL
#  Fly: w1118 Control vs miRNA Mutant
#  Fig: 2

# Rationale
#  To determine whether miRNA mutant line changes FHV-induced mortality compared to paired control.

# Approach
#  Test whether data meets proportional hazards assumption
#  If yes, proceed with cox model as established in (10.1093/evolut/qpae101 and 10.1016/j.virol.2025.110759 )
#  If no, proceed with alternative, RMST (restricted mean survival time, use outlined in 10.3348/kjr.2022.0061 and 10.1186/1471-2288-13-152)

# Preparation
library(survival)
library(survRM2)
library(broom)
library(ggplot2)

# Importing Data
#   Define path to FHV data csvs 
folder_path <- "FHV"
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)


# Test Proportional Hazards Assumption
ph_results <- list()

for (file in file_list) {
  
  # Load data
  data <- read.csv(file, header = TRUE)
  
  # Subset
  subset_data <- subset(data, Virus == 1 & Mutant %in% c(0, 1))
  if (nrow(subset_data) == 0) next
  
  # Fit Cox model
  cox_model <- coxph(Surv(Days, Survival) ~ Mutant, data = subset_data)
  
  # PH test
  ph_test <- cox.zph(cox_model)
  
  # Extract only Mutant row
  ph_df <- data.frame(
    file = basename(file),
    chisq = ph_test$table["Mutant", "chisq"],
    p_value = ph_test$table["Mutant", "p"],
    PH_violation = ph_test$table["Mutant", "p"] < 0.05,
    stringsAsFactors = FALSE
  )
  
  ph_results[[length(ph_results) + 1]] <- ph_df
}

final_ph_results <- do.call(rbind, ph_results)
print(final_ph_results)

# Interpretation: 10/35 Datasets violate the proportional hazards assumption
#  Therefore, cannot proceed with cox model analysis
#  Alternatively, use RMST as an alternative

# Analyse Data using RMST

# Loop over each file
all_results <- list()

for (file in file_list) {
  
  # Load data
  data <- read.csv(file, header = TRUE)
  
  # Subset as required
  subset_data <- subset(data, Virus == 1 & Mutant %in% c(0, 1))
  if (nrow(subset_data) == 0) next
  subset_data$group <- subset_data$Mutant
  
  # Choose tau (minimum of max follow-up time between control/mutant)
  tau <- min(max(subset_data$Days[subset_data$Survival == 1]), 25)
  
  # Run RMST analysis
  rmst_result <- rmst2(
    time   = subset_data$Days,
    status = subset_data$Survival,
    arm    = subset_data$group,
    tau    = tau
  )
  
  # Extract difference row (RMST arm1 - arm0)
  diff_row <- rmst_result$unadjusted.result["RMST (arm=1)-(arm=0)", ]
  
  # Collect results into data frame
  res_df <- data.frame(
    file = basename(file),
    tau = tau,
    RMST_diff = as.numeric(diff_row["Est."]),
    RMST_diff_lower_CI = as.numeric(diff_row["lower .95"]),
    RMST_diff_upper_CI = as.numeric(diff_row["upper .95"]),
    p_value = as.numeric(diff_row["p"]),
    stringsAsFactors = FALSE
  )
  all_results[[length(all_results) + 1]] <- res_df
}
final_results <- do.call(rbind, all_results)
print(final_results)

# Export Results
write.csv(final_results, "rmst_results FHV.csv", row.names = FALSE)

# Next: Plot RMST results for Fig 2
final_results$miR <- sub("\\.csv$", "", final_results$file)
final_results$significance <- ifelse(final_results$p_value <= 0.05,
                                     "Significant (p ≤ 0.05)", 
                                     "Not Significant (p > 0.05)")

final_results$color <- ifelse(final_results$p_value <= 0.05, "#1f78b4", "#e31a1c")

# Reorder by RMST diff. value
final_results$miR <- factor(final_results$miR, 
                            levels = final_results$miR[order(final_results$RMST_diff)])


# Plot Fig. 2
ggplot(final_results, aes(x = RMST_diff, y = miR, color = significance)) +
  geom_errorbar(
    aes(xmin = RMST_diff_lower_CI, xmax = RMST_diff_upper_CI),
    height = 0.25,              
    linewidth = 1.2             
  ) +
  geom_point(size = 3.8, stroke = 1) + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40", linewidth = 1.1) +
  
  scale_color_manual(values = c(
    "Significant (p ≤ 0.05)" = "#1f78b4",  # blue
    "Not Significant (p > 0.05)" = "#e31a1c"  # red
  )) +
  
  scale_x_continuous(
    limits = c(-6, 6),
    breaks = seq(-6, 6, by = 1)
  ) +
  
  labs(
    title = "",
    x = "Restricted Mean Survival Time Difference (FHV)",
    y = "miRNA mutant",
    color = "Significance"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 14, face = "bold"), 
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(face = "bold", size = 14)
  )
