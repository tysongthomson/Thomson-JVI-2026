# Experimental details
#  Virus: PBS only (Mock; collected during DCV replicates)
#  Titre: 1x Phosphate Buffered Saline
#  Fly: w1118 Control vs miRNA Mutant
#  Table S5

# Rationale
#  To determine whether miRNA mutant line changes injection-induced mortality compared to paired control.

# Approach
#  As with DCV, utilise RMST as outlined in 1-2.
#  Notably: Due to low mortality in mock-injected populations, do not include BiolRep as a covariate.
#  Notably: Use same tau as with DCV Analysis, as to check survival difference during experiment timeline..

# References
# 1. Han K, Jung I. 2022. Restricted mean survival time for survival analysis: A quick guide for clinical researchers. Korean J Radiol 23:495.
# 2. Royston P, Parmar MK. 2013. Restricted mean survival time: an alternative to the hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. BMC Med Res Methodol 13:152.

# Preparation
library(survival)
library(survRM2)
library(broom)

# Folder with miRNA CSV files
folder_path <- "DCV"

# Load tau table
tau_table <- read.csv(
  "Table S1. RMST_DCV.csv",
  stringsAsFactors = FALSE
)

# Clean miRNA names
tau_table$miRNAmutant <- trimws(
  tau_table$miRNAmutant
)

# Create tau lookup
tau_lookup <- setNames(
  tau_table$tau,
  tau_table$miRNAmutant
)

# Get CSV files
file_list <- list.files(
  folder_path,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Store results
all_results <- list()

# Loop through files
for (file in file_list) {
  
  # miRNA name
  file_name <- trimws(
    sub("\\.csv$", "", basename(file))
  )
  
  cat("Processing:", file_name, "\n")
  
  # Load data
  data <- read.csv(file)
  
  # Tau for this miRNA
  tau_input <- tau_lookup[[file_name]]
  
  # Skip if no tau
  if (is.null(tau_input) || is.na(tau_input)) next
  
  # Mock data only
  subset_data <- subset(
    data,
    Virus == 0 & Mutant %in% c(0, 1)
  )
  
  # Remove missing values
  subset_data <- subset_data[
    complete.cases(
      subset_data[, c(
        "Days",
        "Survival",
        "Mutant"
      )]
    ),
  ]
  
  # Check mortality in both groups
  event_check <- tapply(
    subset_data$Survival,
    subset_data$Mutant,
    sum
  )
  
  if (
    length(event_check) < 2 ||
    any(event_check == 0)
  ) {
    
    all_results[[length(all_results) + 1]] <- data.frame(
      miRNAmutant = file_name,
      tau = tau_input,
      RMST_diff = NA,
      RMST_diff_lower_CI = NA,
      RMST_diff_upper_CI = NA,
      CI_pm = NA,
      p_value = NA,
      Note =
        "RMST N/A (no mortality for one or both fly lines)"
    ) # Some of the measured fly lines have no mortality recorded for mock-injected. Causing RMST failure.
    
    next
  }
  
  # Safe tau
  arm_max <- tapply(
    subset_data$Days,
    subset_data$Mutant,
    max
  )
  
  tau <- min(
    tau_input,
    min(arm_max)
  )
  
  # RMST analysis
  rmst_result <- tryCatch({
    
    rmst2(
      time = subset_data$Days,
      status = subset_data$Survival,
      arm = subset_data$Mutant,
      tau = tau
    )
    
  }, error = function(e) {
    
    return(NULL)
  })
  
  # Failed RMST
  if (is.null(rmst_result)) {
    
    all_results[[length(all_results) + 1]] <- data.frame(
      miRNAmutant = file_name,
      tau = tau,
      RMST_diff = NA,
      RMST_diff_lower_CI = NA,
      RMST_diff_upper_CI = NA,
      CI_pm = NA,
      p_value = NA,
      Note = "RMST not possible"
    ) # This could occur when there is not enough data.
    
    next
  }
  
  # Extract RMST difference
  diff_table <- rmst_result$unadjusted.result
  
  res <- diff_table["RMST (arm=1)-(arm=0)", ]
  
  ci_pm <- (
    as.numeric(res["upper .95"]) -
      as.numeric(res["lower .95"])
  ) / 2
  
  # Save successful result
  all_results[[length(all_results) + 1]] <- data.frame(
    miRNAmutant = file_name,
    tau = tau,
    RMST_diff = as.numeric(res["Est."]),
    RMST_diff_lower_CI =
      as.numeric(res["lower .95"]),
    RMST_diff_upper_CI =
      as.numeric(res["upper .95"]),
    CI_pm = ci_pm,
    p_value = as.numeric(res["p"]),
    Note = ""
  )
}

# Combine results
final_results <- do.call(
  rbind,
  all_results
)

# Sort by miRNA name
final_results <- final_results[
  order(final_results$miRNAmutant),
]

# Save output
write.csv(
  final_results,
  "Table S6. RMST_DCV Mock-only.csv",
  row.names = FALSE
)

# Print results
print(final_results)
