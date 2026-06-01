# Experimental details
#  Virus: Drosophila C virus (DCV, isolate EB) [RefSeq: NC_001834.1]
#  Titre: 10E7 IU/mL
#  Fly: w1118 Control vs miRNA Mutant
#  Fig. 1, Table S1

# Rationale
#  To determine whether miRNA mutant line changes DCV-induced mortality compared to paired control.

# Approach
#  Test whether data meets proportional hazards assumption
#  If yes, proceed with cox model as established in (1-2)
#  If no, proceed with alternative, RMST (restricted mean survival time, use outlined in 3-4)

# References
# 1. Mendel BM, Asselin AK, Johnson KN, McGuigan K. 2024. Effects of spontaneous mutations on survival and reproduction of Drosophila serrata infected with Drosophila C virus. Evolution 78:1661–1672.
# 2. Therneau T. 2024. coxme: Mixed effects cox models. CRAN Repos.
# 3. Han K, Jung I. 2022. Restricted mean survival time for survival analysis: A quick guide for clinical researchers. Korean J Radiol 23:495.
# 4. Royston P, Parmar MK. 2013. Restricted mean survival time: an alternative to the hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. BMC Med Res Methodol 13:152.

# Preparation
library(survival)
library(survRM2)
library(broom)
library(ggplot2)

# Importing Data
#   Define path to DCV data csvs 
folder_path <- "DCV"
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

# Interpretation: 11/34 Datasets violate the proportional hazards assumption
#  Therefore, cannot proceed with cox model analysis
#  As such, use RMST as an alternative

# Analyse Data using RMST

# Define folder containing miRNA CSV files
folder_path <- "DCV"

# Get list of all CSV files in folder
file_list <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

# Prepare list to store results from each file
all_results <- list()

# Loop over each miRNA dataset
for (file in file_list) {
  
  cat("Processing:", basename(file), "\n")
  
  # Load data
  data <- read.csv(file)
  
  # Subset to infection condition and relevant groups
  subset_data <- subset(data, Virus == 1 & Mutant %in% c(0,1))
  
  # Remove incomplete cases to prevent covariate mismatch errors (shouldn't be necessary)
  subset_data <- subset_data[
    complete.cases(subset_data[, c("Days", "Survival", "Mutant", "BiolRep")]),
  ]
  
  # Skip if no data after filtering
  if (nrow(subset_data) == 0) next
  
  # Define treatment group variable
  subset_data$group <- subset_data$Mutant
  
  # Convert biological replicate to factor
  subset_data$BiolRep <- factor(subset_data$BiolRep)
  
  # Skip if only one biological replicate exists
  if (length(unique(subset_data$BiolRep)) < 2) next
  
  # Create covariate matrix for adjustment
  covars <- model.matrix(~ BiolRep, data = subset_data)[, -1, drop = FALSE]
  
  # Define tau (maximum observed event time)
  tau <- max(subset_data$Days[subset_data$Survival == 1])
  
  # Run RMST analysis with biological replicate adjustment
  rmst_result <- rmst2(
    time = subset_data$Days,
    status = subset_data$Survival,
    arm = subset_data$group,
    tau = tau,
    covariates = covars
  )
  
  # Extract adjusted RMST results
  diff_table <- rmst_result$adjusted.result
  
  # Identify correct RMST comparison row
  row_name <- grep("arm=1", rownames(diff_table), value = TRUE)[1]
  
  res <- diff_table[row_name, ]
  
  # Compute CI half-width
  ci_pm <- (as.numeric(res["upper .95"]) - as.numeric(res["lower .95"])) / 2
  
  # Store results (remove .csv from filename)
  res_df <- data.frame(
    miRNAmutant = sub("\\.csv$", "", basename(file)),
    tau = tau,
    RMST_diff = as.numeric(res["Est."]),
    RMST_diff_lower_CI = as.numeric(res["lower .95"]),
    RMST_diff_upper_CI = as.numeric(res["upper .95"]),
    CI_pm = ci_pm,
    p_value = as.numeric(res["p"])
  )
  
  # Append to results list
  all_results[[length(all_results) + 1]] <- res_df
}

# Combine all miRNA results into one table
final_results <- do.call(rbind, all_results)

# Save output table
write.csv(final_results, "Table S1. RMST_DCV.csv", row.names = FALSE)

# Print final results
print(final_results)

# Next: Plot RMST results for Fig 1
final_results$miR <- sub("\\.csv$", "", final_results$miRNAmutant)
final_results$significance <- ifelse(final_results$p_value <= 0.05,
                                     "Significant (p ≤ 0.05)", 
                                     "Not Significant (p > 0.05)")

final_results$color <- ifelse(final_results$p_value <= 0.05, "#1f78b4", "#e31a1c")

# Reorder by RMST diff. value
final_results$miR <- factor(final_results$miR, 
                            levels = final_results$miR[order(final_results$RMST_diff)])


# Plot Fig. 1
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
    x = "Mean Survival Time Difference (DCV)",
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

# Save Fig. 1
ggsave(
  filename = "Fig. 1.tiff",
  plot = last_plot(),
  device = "tiff",
  units = "px",
  width = 5000, 
  height = 3000,
  dpi = 300,
  compression = "lzw"
)
