# RScript for Generating the Survival Curves utilised for Figure S1
# This script utilises the DCV folder data and rmst_results DCV.csv provided in Figure 1 folder.

# Preparation
library(survival)
library(survminer)
library(ggplot2)

# Importing Data/Defining Paths
data_dir <- "DCV"
output_dir <- file.path(data_dir, "survival_plots")
dir.create(output_dir, showWarnings = FALSE)

# Load RMST results
rmst <- read.csv(file.path("rmst_results DCV.csv"), stringsAsFactors = FALSE)

# Get all CSV files in the data directory
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

# Generate Survival Curves (Fig. S1)
for (file in files) {
  
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Only keep Virus == 1 (censor PBS control)
  df <- df[df$Virus == 1, ]
  if (nrow(df) == 0) next
  
  # Extract miRNA name from filename
  file_base <- tools::file_path_sans_ext(basename(file))
  miRNA_name <- strsplit(file_base, "_")[[1]][1]
  
  # Get RMST, p-value, tau results for each miRNA
  rmst_row <- rmst[rmst$file == paste0(miRNA_name, ".csv"), ]
  
  rmst_text <- NULL
  if (nrow(rmst_row) == 1) {
    
    rmst_diff <- round(rmst_row$RMST_diff, 2)
    
    if (rmst_row$p_value < 0.05) {
      p_text <- "p <= 0.05"
    } else {
      p_text <- paste0("p = ", signif(rmst_row$p_value, 2))
    }
    
    tau_text <- paste0("tau = ", round(rmst_row$tau, 2))
    
    rmst_text <- paste0(
      "RMST diff = ", rmst_diff, "\n",
      p_text, "\n",
      tau_text
    )
  }
  
  df$Mutant <- factor(
    df$Mutant,
    levels = c(0, 1),
    labels = c("Control", paste0(miRNA_name, "-KO"))
  )
  
  # Survival object and fit
  surv_obj <- Surv(time = df$Days, event = df$Survival)
  fit <- survfit(surv_obj ~ Mutant, data = df)
  
  xmax <- max(df$Days[df$Survival == 1], na.rm = TRUE)
  if (!is.finite(xmax)) xmax <- max(df$Days)
  
  # Create survival plot
  p <- ggsurvplot(
    fit,
    data = df,
    conf.int = FALSE,
    censor = TRUE,
    risk.table = FALSE,
    xlim = c(1, xmax),
    break.time.by = 1,
    xlab = "Days",
    ylab = "Proportional Survival (DCV)",
    legend.title = "Condition",
    legend.labs = c("Control", paste0(miRNA_name, "-KO")),
    ggtheme = theme_bw() +
      theme(
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 12, color = "black"),
        axis.text.y  = element_text(size = 12, color = "black"),
        legend.text  = element_text(size = 12)
      ),
    palette = c("blue", "red")
  )
  
  # Set axes scales and legend position
  p$plot <- p$plot +
    scale_x_continuous(
      breaks = seq(0, xmax, by = 1),
      limits = c(0, xmax),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    theme(legend.position = "right")
  
  # Add RMST annotation if available
  if (!is.null(rmst_text)) {
    # Dynamic y-position: 10% below max survival
    y_pos <- 0.1 + 0.05
    p$plot <- p$plot +
      annotate(
        "text",
        x = 2,          # x-position (days)
        y = y_pos,      # y-position
        label = rmst_text,
        hjust = 0,
        size = 4
      )
  }
  
  # Save plot
  out_file <- file.path(
    output_dir,
    paste0(file_base, "_DCV_survival_curve.pdf")
  )
  
  ggsave(
    filename = out_file,
    plot = p$plot,
    width = 10,
    height = 6,
    dpi = 300,
    device = "pdf"
  )

}

# Individual Survival Curves were combined using external software
