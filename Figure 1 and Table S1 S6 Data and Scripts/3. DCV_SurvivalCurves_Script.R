# Experimental details
#  Virus: Drosophila C virus (DCV, isolate EB) [RefSeq: NC_001834.1]
#  Titre: 10E7 IU/mL
#  Fly: w1118 Control vs miRNA Mutant
#  Fig. S1

# Rationale
#  Generate survival curves for screened lines, allowing for visualisation of miRNA mutant vs control comparison.
#  Additionally, include mock-injected flies.

# Preparation
library(survival)
library(survminer)
library(ggplot2)

# Importing Data / Defining Paths
data_dir <- "DCV"

output_dir <- file.path("Fig. S1 survival_plots")
dir.create(output_dir, showWarnings = FALSE)

# Load DCV RMST results
rmst <- read.csv(
  "Table S1. RMST_DCV.csv",
  stringsAsFactors = FALSE
)

# Load Mock-only RMST results
rmst_mock <- read.csv(
  "Table S6. RMST_DCV Mock-only.csv",
  stringsAsFactors = FALSE
)

# Get all CSV files
files <- list.files(
  data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Generate Survival Curves
for (file in files) {
  
  # Load dataset
  df <- read.csv(file, stringsAsFactors = FALSE)
  if (nrow(df) == 0) next
  
  # Extract miRNA name
  file_base <- tools::file_path_sans_ext(basename(file))
  miRNA_name <- strsplit(file_base, "_")[[1]][1]
  
  # RMST lookup
  rmst_row <- rmst[rmst$miRNAmutant == miRNA_name, ]
  rmst_mock_row <- rmst_mock[rmst_mock$miRNAmutant == miRNA_name, ]
  
  rmst_text <- NULL
  rmst_mock_text <- NULL
  
  # DCV RMST (assumed clean)
  if (nrow(rmst_row) == 1) {
    
    rmst_diff <- round(rmst_row$RMST_diff, 2)
    
    if (rmst_row$p_value < 0.05) {
      p_text <- "p <= 0.05"
    } else {
      p_text <- paste0("p = ", signif(rmst_row$p_value, 2))
    }
    
    tau_text <- paste0("tau = ", round(rmst_row$tau, 2))
    
    rmst_text <- paste(
      "DCV",
      paste0("RMST diff = ", rmst_diff),
      p_text,
      tau_text,
      sep = "\n"
    )
    
    xmax <- ceiling(rmst_row$tau + 1)
    
  } else {
    
    xmax <- max(df$Days, na.rm = TRUE)
  }
  
  # Mock RMST
  if (nrow(rmst_mock_row) == 1 &&
      !all(is.na(rmst_mock_row$RMST_diff)) &&
      !all(is.na(rmst_mock_row$p_value)) &&
      !all(is.na(rmst_mock_row$tau))) {
    
    mock_diff <- round(rmst_mock_row$RMST_diff, 2)
    
    if (rmst_mock_row$p_value < 0.05) {
      mock_p <- "p <= 0.05"
    } else {
      mock_p <- paste0("p = ", signif(rmst_mock_row$p_value, 2))
    }
    
    rmst_mock_text <- paste(
      "Mock",
      paste0("RMST diff = ", mock_diff),
      mock_p,
      sep = "\n"
    )
    
  } else {
    
    rmst_mock_text <- paste(
      "Mock",
      "RMST diff = NA",
      "p = NA",
      sep = "\n"
    )
  }
  
  # Create combined group labels
  df$Group <- with(
    df,
    paste(
      ifelse(Mutant == 0,
             "Control",
             paste0(miRNA_name, "-KO")),
      ifelse(Virus == 0,
             "Mock",
             "DCV"),
      sep = "-"
    )
  )
  
  df$Group <- factor(
    df$Group,
    levels = c(
      "Control-Mock",
      paste0(miRNA_name, "-KO-Mock"),
      "Control-DCV",
      paste0(miRNA_name, "-KO-DCV")
    )
  )
  
  # Survival analysis
  surv_obj <- Surv(time = df$Days, event = df$Survival)
  fit <- survfit(surv_obj ~ Group, data = df)
  
  # Plot survival curves
  p <- ggsurvplot(
    
    fit,
    data = df,
    
    conf.int = FALSE,
    censor = TRUE,
    risk.table = FALSE,
    
    break.time.by = 1,
    
    xlab = "Days",
    ylab = "Proportional Survival",
    
    legend.title = "Condition",
    legend.labs = levels(df$Group),
    
    palette = c("grey40", "grey70", "blue", "red"),
    linetype = c("dashed", "longdash", "solid", "twodash"),
    linewidth = 1.2,
    
    ggtheme = theme_bw() +
      theme(
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 12, color = "black"),
        axis.text.y  = element_text(size = 12, color = "black"),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        legend.position = "right"
      )
  )
  
  # Axis formatting
  p$plot <- p$plot +
    
    scale_x_continuous(
      breaks = seq(0, 25, by = 1),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, 1.01),
      breaks = seq(0, 1, by = 0.1),
      expand = c(0, 0)
    ) +
    
    coord_cartesian(xlim = c(1, xmax))
  
  # Add RMST annotations (Mock above DCV)
  if (!is.null(rmst_mock_text)) {
    p$plot <- p$plot +
      annotate(
        "text",
        x = 1.5,
        y = 0.40,
        label = rmst_mock_text,
        hjust = 0,
        size = 4,
        fontface = "bold"
      )
  }
  
  if (!is.null(rmst_text)) {
    p$plot <- p$plot +
      annotate(
        "text",
        x = 1.5,
        y = 0.15,
        label = rmst_text,
        hjust = 0,
        size = 4,
        fontface = "bold"
      )
  }
  
  # Save plot
  out_file <- file.path(
    output_dir,
    paste0(file_base, "_survival_curve.pdf")
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

