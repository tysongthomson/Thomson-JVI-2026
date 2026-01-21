# Analysis details

# Rationale/Questions
#  We want to use our results to determine whether there is a correlation between conservation and virus impact.
#  To address this rationale, we had the following questions:
#  1. Are miRNAs which impact virus infection are more likely to be conserved (than no-impact)?
#  2. Are conserved miRNAs more likely to target more/multiple viruses (than non-conserved)?

# Approach
#  To address Question 1.
#  a) Fisher's exact test: Comparing Conservation and Impact on virus (binary).
#  b) Linear Model: Comparing Conservation and Impact on virus (continuous).
#  c) Wilcoxon rank-sum test: Comparing Conservation and Impact on virus (continuous).
#  d) Fisher's exact test: Comparing Conservation and Impact on virus (categorical).
#  e) Linear Model: Comparing Conservation and Impact on virus (continuous, absolute).
#
#  To address Question 2.
#  a) Ordinal Logistic Model: Comparing Conservation and Number of Viruses impacted (categorical).
#  b) Wilcoxon rank-sum test: Comparing Conservation and Number of Viruses impacted (categorical).

# Preparation
library(Biostrings)
library(dplyr)
library(tidyr)
library(MASS)

# Importing Data
mirnas <- readRNAStringSet("allmirnas.fa") # all miRNA seqs from miRBase

# Step 0: Extract mature miRNA sequences.
df <- data.frame(
  header = names(mirnas),
  sequence = as.character(mirnas),
  stringsAsFactors = FALSE
)

df <- df[grepl("^dme-", df$header), ]

# All screened miRNAs
target_mirnas <- c(
  "miR-10","miR-1007","miR-1014","miR-11","miR-124","miR-133","miR-137",
  "miR-219","miR-263b","miR-274","miR-276b","miR-278","miR-282","miR-283",
  "miR-285","miR-2b-1","miR-2b","miR-304","miR-31a","miR-31b","miR-375","miR-87",
  "miR-92a","miR-932","miR-957","miR-958","miR-965","miR-966","miR-967",
  "miR-971","miR-986","miR-988","miR-989","miR-990","miR-999","miR-9c"
)

pattern <- paste0("dme-(", paste(target_mirnas, collapse = "|"), ")-[35]p")

df_target <- df[grepl(pattern, df$header), ]

df_target$mirna_id <- sub("^([^ ]+).*", "\\1", df_target$header)
df_target$mimat_id <- sub("^[^ ]+ ([^ ]+).*", "\\1", df_target$header)
df_target$arm <- ifelse(grepl("-5p", df_target$mirna_id), "5p", "3p")
df_target$mirna_name <- sub("^dme-(miR-[^-]+).*", "\\1", df_target$mirna_id)

final_table <- df_target[, c(
  "mirna_name","mirna_id","arm","mimat_id","sequence"
)]

print(final_table)

rna_out <- RNAStringSet(final_table$sequence)
names(rna_out) <- final_table$mirna_id
# Save extracted miRNA sequences
writeXStringSet(rna_out, "dme_selected_mature_miRNAs.fa")

# Importing Data
dmemirs <- readRNAStringSet("dme_selected_mature_miRNAs.fa")
allmir <- mirnas

# Step 0.5: Establish Conservation levels between Dme and Dsi, Aae and Aga.
#  Note: To be identified as conserved, the matching miRNA required a 100% seed region match and ≥80% sequence similarity.
#        Additionally, for each miRNA gene, the most conserved mature miRNA arm was selected.


df_all <- data.frame(
  name = names(allmir),
  seq = as.character(allmir),
  stringsAsFactors = FALSE
)

species <- c("dsi", "aae", "aga") # D. simulans, A. aegypti and A. gambiae.

get_seed <- function(seq) substr(seq, 2, 8) # Specify seed region of miRNAs

pct_identity <- function(a, b) {
  len <- min(nchar(a), nchar(b))
  sum(strsplit(substr(a,1,len),"")[[1]] ==
        strsplit(substr(b,1,len),"")[[1]]) / len * 100
}

results <- list()

# Loops through and queries all Dme miRNAs, identifying the highest conservation miRNA for each of Dsi, Aae and Aga
for (qname in names(dmemirs)) {
  qseq <- as.character(dmemirs[[qname]])
  qseed <- get_seed(qseq)
  
  for (sp in species) {
    candidates <- df_all[grepl(paste0("^", sp, "-"), df_all$name), ]
    candidates$seed <- substr(candidates$seq, 2, 8)
    candidates <- candidates[candidates$seed == qseed, ]
    
    if (nrow(candidates) == 0) {
      results[[length(results)+1]] <- data.frame(
        dmemirs = qname,
        species = sp,
        best_match = NA,
        pct_identity = NA,
        stringsAsFactors = FALSE
      )
    } else {
      candidates$pct <- sapply(candidates$seq, pct_identity, b = qseq)
      best <- candidates[which.max(candidates$pct), ]
      
      results[[length(results)+1]] <- data.frame(
        dmemirs = qname,
        species = sp,
        best_match = best$name,
        pct_identity = best$pct,
        stringsAsFactors = FALSE
      )
    }
  }
}

final <- do.call(rbind, results)
print(final)
write.csv(final, "miRNA_conservation.csv", row.names = FALSE)


# Importing Data
cons <- final

cons$conserved <- ifelse(!is.na(cons$pct_identity) & cons$pct_identity >= 80, TRUE, FALSE)

cons <- cons %>%
  mutate(
    miRNA_clean = gsub("^dme-", "", dmemirs),
    miRNA_base  = gsub("-(3p|5p)$", "", miRNA_clean)
  )

# Collated results is the RMST/T-Test results from each Table S1-3, compiled into a single csv.
results <- read.csv("collated_results.csv", stringsAsFactors = FALSE)

# Step 0.75: Create a collated file for all Virus impacts, and conservation levels.
#   Collapse conservation per miRNA (3p/5p combined)
#   If only one of -3p and -5p are considered conserved, the miRNA is still considered conserved.
cons_summary <- cons %>%
  group_by(miRNA_base, species) %>%
  summarise(
    conserved = any(conserved),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = species,
    values_from = conserved,
    names_prefix = "conserved_",
    values_fill = FALSE
  )

# Add virus impact columns
results <- results %>%
  mutate(
    DCV_impact  = DCV_p_value < 0.05,
    FHV_impact  = FHV_p_value < 0.05,
    SINV_impact = ifelse(is.na(SINV_p_value), TRUE, SINV_p_value < 0.05)
  )

# Merge conservation with results
final_table <- results %>%
  left_join(cons_summary, by = c("miRNA" = "miRNA_base"))

print(final_table)
write.csv(final_table, "miRNA_data_collated.csv", row.names = FALSE)

# Step 1: Address Question 1

# Step 1a: Fisher's exact test: Comparing Conservation and Impact on virus (binary).
#  Data Structure:
#  Conservation: TRUE/FALSE binary for each of Dsi, Aae and Aga
#  Virus Impact: TRUE/FALSE binary for each of DCV, FHV and SINV

# Null Hypothesis: There is no correlation between whether miRNA impacts virus and conservation status.
# Note: This analysis treats Impact as binary (doesn't distingush between increase/decrease)

# Function to run Fisher's exact test
run_fisher <- function(impact_col, conserved_col, data) {
  tab <- table(
    Impact = data[[impact_col]],
    Conserved = data[[conserved_col]]
  )
  
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0]
  
  if(all(dim(tab) == c(2,2))) {
    test <- fisher.test(tab)
    list(odds_ratio = test$estimate, p_value = test$p.value)
  } else {
    list(odds_ratio = NA, p_value = NA)
  }
}

viruses <- c("DCV", "FHV", "SINV")
conservation <- c("dsi", "aae", "aga")

results_list <- list()
for(v in viruses){
  for(c in conservation){
    impact_col <- paste0(v, "_impact")
    conserved_col <- paste0("conserved_", c)
    results_list[[paste(v, c, sep = "_")]] <- run_fisher(impact_col, conserved_col, final_table)
  }
}

# Convert list to data frame
results_df <- do.call(rbind, lapply(names(results_list), function(x){
  data.frame(
    test = x,
    odds_ratio = results_list[[x]]$odds_ratio,
    p_value = results_list[[x]]$p_value
  )
}))

results_df$odds_ratio <- as.numeric(results_df$odds_ratio)
results_df$p_value <- as.numeric(results_df$p_value)

print(results_df)
write.csv(results_df, "Question 1a.csv", row.names = FALSE)

# Interpretation
#  P-value is above 0.05 for all, cannot refute Null hypothesis.
#  Therefore, there is no evidence that conservation and miRNA impact on virus is correlated.
#  This was consistent for all virus and species combinations.

# Step 1b: Linear Model: Comparing Conservation and Impact on virus (continuous).
#  Data Structure:
#  Conservation: TRUE/FALSE binary for each of Dsi, Aae and Aga
#  Virus Impact: Continuous RMST diff/Log2FC for each of DCV, FHV and SINV

# Null Hypothesis: Conservation status has no impact on magnitude of virus impact.
# Note: This analysis treats Impact as continuous (doesn't consider p-value of virus impact)

# List of viruses and conservation types
viruses <- c("DCV", "FHV", "SINV")
conservation <- c("dsi", "aae", "aga")

# Initialize results data frame
results_lm <- data.frame(
  test = character(),
  estimate = numeric(),
  p_value = numeric(), 
  stringsAsFactors = FALSE
)

# Loop through each virus and conservation type
for(v in viruses){
  for(c in conservation){
    
    # Determine which column holds the viral impact
    # RMST_diff for DCV/FHV, log2 fold-change for SINV
    impact_col <- if(v %in% c("DCV","FHV")) paste0(v, "_RMST_diff") else "SINV_log2_fold_change"
    conserved_col <- paste0("conserved_", c)
    
    # Subset to rows where the impact is not NA
    df <- final_table[!is.na(final_table[[impact_col]]), ]
    
    # Ensure the conservation column is treated as a factor
    df[[conserved_col]] <- as.factor(df[[conserved_col]])
    
    # Fit linear model: impact ~ conservation
    lm_fit <- lm(df[[impact_col]] ~ df[[conserved_col]])
    
    coef_name <- grep("TRUE", names(coef(lm_fit)), value = TRUE)
    
    if(length(coef_name) == 1){
      # Extract estimate and p-value for conserved = TRUE
      est <- coef(lm_fit)[coef_name]
      pval <- summary(lm_fit)$coefficients[coef_name, "Pr(>|t|)"]
    } else {
      est <- NA
      pval <- NA
    }
    
    # Add results to final table
    results_lm <- rbind(
      results_lm,
      data.frame(
        test = paste(v, c, sep = "_"),
        estimate = est,
        p_value = pval,
        stringsAsFactors = FALSE
      ),
      row.names = NULL
    )
  }
}

# Print and save results
print(results_lm)
write.csv(results_lm, "Question 1b.csv", row.names = FALSE)

# Interpretation
#  P-value is above 0.05 for all, cannot refute Null hypothesis.
#  Therefore, there is no evidence that conservation and magnitude of miRNA impact on virus is correlated.
#  This was consistent for all virus and species combinations.

# Step 1c: Wilcoxon rank-sum test: Comparing Conservation and Impact on virus distribution (continous).
#  Data Structure:
#  Conservation: TRUE/FALSE binary for each of Dsi, Aae and Aga
#  Virus Impact: Continuous RMST diff/Log2FC for each of DCV, FHV and SINV

# Null Hypothesis: Conservation status does not impact distribution of virus impacts.
# Note: This analysis treats Impact as continuous (doesn't consider p-value of virus impact)

viruses <- c("DCV", "FHV", "SINV")
conservation <- c("dsi", "aae", "aga")

# Initialize results data frame
results_wilcox <- data.frame(
  test = character(),       # virus_conservation combination
  W_statistic = numeric(),  # Wilcoxon W statistic
  p_value = numeric(),      # p-value
  stringsAsFactors = FALSE
)

# Loop through each virus and conservation type
for(v in viruses){
  for(c in conservation){
    
    # Determine which column holds the viral impact
    # RMST_diff for DCV/FHV, log2 fold-change for SINV
    impact_col <- if(v %in% c("DCV","FHV")) paste0(v, "_RMST_diff") else "SINV_log2_fold_change"
    conserved_col <- paste0("conserved_", c)
    
    # Subset rows with non-NA impact and non-NA conservation
    df <- final_table[!is.na(final_table[[impact_col]]) & !is.na(final_table[[conserved_col]]), ]
    
    # Ensure conservation is a factor (TRUE/FALSE)
    df[[conserved_col]] <- as.factor(df[[conserved_col]])
    
    # Wilcoxon rank-sum test: impact ~ conservation
    wilcox_res <- wilcox.test(df[[impact_col]] ~ df[[conserved_col]], exact = FALSE)
    
    results_wilcox <- rbind(
      results_wilcox,
      data.frame(
        test = paste(v, c, sep = "_"),
        W_statistic = wilcox_res$statistic,
        p_value = wilcox_res$p.value,
        stringsAsFactors = FALSE
      ),
      row.names = NULL
    )
  }
}

print(results_wilcox)
write.csv(results_wilcox, "Question 1c.csv", row.names = FALSE)

# Interpretation
#  P-value is above 0.05 for all, cannot refute Null hypothesis.
#  Therefore, there is no evidence that conservation status changes distribution of virus impacts.
#  This was consistent for all virus and species combinations.

# Step 1d: Fisher's exact test: Comparing Conservation and Impact on virus (categorical).
#  Data Structure:
#  Conservation: TRUE/FALSE binary for each of Dsi, Aae and Aga
#  Virus Impact: POS/NEG/FALSE categorical for each of DCV, FHV and SINV

# Null Hypothesis: There is no correlation between whether miRNA impacts virus and conservation status.
# Note: This analysis distinguished between POS/NEG impacts on virus.

# Define viruses and conservation types
viruses <- c("DCV", "FHV", "SINV")
conservation <- c("dsi", "aae", "aga")

# Initialize results table
results_fisher3x2 <- data.frame(
  test = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(v in viruses){
  # Determine which column holds effect size
  effect_col <- if(v %in% c("DCV","FHV")) paste0(v, "_RMST_diff") else "SINV_log2_fold_change"
  impact_col <- paste0(v, "_impact")
  
  for(c in conservation){
    conserved_col <- paste0("conserved_", c)
    
    # Create impact category column
    impact_cat <- rep("FALSE", nrow(final_table))
    
    # Assign POS/NEG for TRUE impacts based on effect size
    idx_true <- which(final_table[[impact_col]] == TRUE)
    impact_cat[idx_true] <- ifelse(final_table[[effect_col]][idx_true] > 0, "POS", "NEG")
    
    # Convert to factor
    impact_cat <- factor(impact_cat, levels = c("NEG", "FALSE", "POS"))
    conserved_fac <- factor(final_table[[conserved_col]], levels = c(FALSE, TRUE))
    
    # Create 3x2 table
    tab <- table(Impact = impact_cat, Conserved = conserved_fac)
    
    # Run Fisher's exact test
    fisher_res <- fisher.test(tab)
    
    # Store results
    results_fisher3x2 <- rbind(
      results_fisher3x2,
      data.frame(
        test = paste(v, c, sep = "_"),
        p_value = fisher_res$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
}

# View and save results
print(results_fisher3x2)
write.csv(results_fisher3x2, "Question 1d.csv", row.names = FALSE)

# Interpretation
#  P-value is above 0.05 for all, cannot refute Null hypothesis.
#  Therefore, there is no evidence that conservation status has an effect on whether miRNA impacts virus.
#  This was consistent for all virus and species combinations.

# Step 1e: Linear Model: Comparing Conservation and Impact on virus (continuous, absolute).
#  Data Structure:
#  Conservation: TRUE/FALSE binary for each of Dsi, Aae and Aga
#  Virus Impact: Continuous RMST diff/Log2FC for each of DCV, FHV and SINV (absolute value)

# Null Hypothesis: Conservation status has no impact on magnitude of virus impact.
# Note: This analysis uses absolute value of impact, does not consider pos/neg.

# List of viruses and conservation types
viruses <- c("DCV", "FHV", "SINV")
conservation <- c("dsi", "aae", "aga")

# Initialize results data frame
results_lm <- data.frame(
  test = character(),
  estimate = numeric(),
  p_value = numeric(), 
  stringsAsFactors = FALSE
)

for(v in viruses){
  for(c in conservation){
    
    impact_col <- if(v %in% c("DCV","FHV")) paste0(v, "_RMST_diff") else "SINV_log2_fold_change"
    conserved_col <- paste0("conserved_", c)
    
    # Subset rows with non-NA impact
    df <- final_table[!is.na(final_table[[impact_col]]), ]
    
    # Absolute value of impact
    df$impact_abs <- abs(df[[impact_col]])
    
    # Ensure conservation is a factor
    df[[conserved_col]] <- as.factor(df[[conserved_col]])
    
    # Correct linear model: specify data = df
    lm_fit <- lm(impact_abs ~ df[[conserved_col]], data = df)
    
    # Extract coefficient for TRUE
    coef_name <- grep("TRUE", names(coef(lm_fit)), value = TRUE)
    
    if(length(coef_name) == 1){
      est <- coef(lm_fit)[coef_name]
      pval <- summary(lm_fit)$coefficients[coef_name, "Pr(>|t|)"]
    } else {
      est <- NA
      pval <- NA
    }
    
    # Add to results
    results_lm <- rbind(
      results_lm,
      data.frame(
        test = paste(v, c, sep = "_"),
        estimate = est,
        p_value = pval,
        stringsAsFactors = FALSE
      ),
      row.names = NULL
    )
  }
}

# View and save
print(results_lm)
write.csv(results_lm, "Question 1e.csv", row.names = FALSE)

# Interpretation
#  P-value is above 0.05 for all, cannot refute Null hypothesis.
#  Therefore, there is no evidence that conservation level and magnitude of miRNA impact on virus is correlated.
#  This was consistent for all virus and species combinations.

# Step 2: Address Question 2

# Step 2a: Ordinal Logistic Model: Comparing Conservation and Number of Viruses impacted (categorical).
#  Data Structure:
#  Conservation: TRUE/FALSE binary for each of Dsi, Aae and Aga
#  Virus Impact: Number of viruses impacted, Categorical (0, 1, 2, 3).

# Null Hypothesis: Conservation status has no effect on number of viruses a miRNA impacts.

# Create a column for number of viruses impacted (0-3)
final_table$viruses_impacted <- rowSums(final_table[, c("DCV_impact", "FHV_impact", "SINV_impact")], na.rm = TRUE)

# Convert to ordered factor (0 < 1 < 2 < 3)
final_table$viruses_impacted <- factor(final_table$viruses_impacted, levels = 0:3, ordered = TRUE)

# Define conservation types
conservation <- c("dsi", "aae", "aga")

# Initialize results table
results_ordinal <- data.frame(
  test = character(),   # conservation type
  estimate = numeric(), # coefficient for conserved = TRUE
  p_value = numeric(),  # Wald test p-value
  stringsAsFactors = FALSE
)

# Loop through conservation types and fit ordinal logistic regression
for(c in conservation){
  
  conserved_col <- paste0("conserved_", c)
  
  # Ensure predictor is factor
  final_table[[conserved_col]] <- as.factor(final_table[[conserved_col]])
  
  # Fit ordinal logistic regression: viruses_impacted ~ conserved
  model <- polr(viruses_impacted ~ final_table[[conserved_col]], data = final_table, Hess = TRUE)
  
  # Extract coefficient name for TRUE level
  coef_name <- grep("TRUE", names(coef(model)), value = TRUE)
  est <- coef(model)[coef_name]
  
  # Calculate p-value using Wald test: z = estimate / SE
  se <- sqrt(diag(vcov(model)))[coef_name]
  z <- est / se
  pval <- 2 * (1 - pnorm(abs(z)))
  
  # Store results
  results_ordinal <- rbind(
    results_ordinal,
    data.frame(
      test = paste("viruses_impacted", c, sep = "_"),
      estimate = est,
      p_value = pval,
      stringsAsFactors = FALSE
    ),
    row.names = NULL
  )
}

# Print and save results
print(results_ordinal)
write.csv(results_ordinal, "Question 2a.csv", row.names = FALSE)

# Interpretation
#  P-value is above 0.05 for all, cannot refute Null hypothesis.
#  Therefore, there is no evidence that conservation status has an effect on number of viruses a miRNA impacts.
#  This was consistent for all species tested.

# Step 2b: Wilcoxon rank-sum test: Comparing Conservation and Number of Viruses impacted (categorical).
#  Data Structure:
#  Conservation: TRUE/FALSE binary for each of Dsi, Aae and Aga
#  Virus Impact: Number of viruses impacted, Categorical (0, 1, 2, 3).

# Null Hypothesis: Conservation status has no effect on distribution of number of viruses impacted.

# Define conservation types
conservation <- c("dsi", "aae", "aga")

# Initialize results table
results_wilcox <- data.frame(
  test = character(),
  W_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each conservation type
for(c in conservation){
  
  conserved_col <- paste0("conserved_", c)
  
  # Subset data for non-NA values of conservation and viruses_impacted
  df <- final_table[!is.na(final_table[[conserved_col]]) & !is.na(final_table$viruses_impacted), ]
  
  # Ensure conservation is factor
  df[[conserved_col]] <- as.factor(df[[conserved_col]])
  
  # Convert ordered factor to numeric for Wilcoxon test
  impact_numeric <- as.numeric(as.character(df$viruses_impacted))
  
  # Perform Wilcoxon rank-sum test: number of viruses impacted ~ conserved
  wilcox_res <- wilcox.test(impact_numeric ~ df[[conserved_col]], exact = FALSE)
  
  # Store results
  results_wilcox <- rbind(
    results_wilcox,
    data.frame(
      test = paste("viruses_impacted", c, sep = "_"),
      W_statistic = wilcox_res$statistic,
      p_value = wilcox_res$p.value,
      stringsAsFactors = FALSE
    ),
    row.names = NULL
  )
}

# Print and save results
print(results_wilcox)
write.csv(results_wilcox, "Question 2b.csv", row.names = FALSE)

# Interpretation
#  P-value is above 0.05 for all, cannot refute Null hypothesis.
#  Therefore, there is no evidence that conservation status has an effect on distribution of number of viruses impacted.
#  This was consistent for all species tested.


