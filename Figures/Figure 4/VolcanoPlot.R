# ============================================================
# Volcano Plot
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

# Load required libraries
library(tidyverse)   # For data manipulation and visualization
library(caret)       # For machine learning functions
library(glmnet)      # For LASSO regression

# List of additional packages needed
packages <- c('reshape2','reshape','plyr','dplyr','tidyr','mice','ggplot2','knitr',
              'lubridate','openair','gridExtra','zoo','car','ggsignif','haven',
              'labelled','cmprsk','riskRegression','foreign','stats','data.table',
              'tsModel','stringr','ggthemes','cowplot','dlnm','scales','splines',
              'mgcv','HH','Epi','FSA','readxl','writexl','compareGroups','DescTools')

# Load all packages
lapply(packages, require, character.only = TRUE)

# Avoid namespace conflicts with other packages
select = dplyr::select
rename = dplyr::rename

# Load training and testing datasets
data1 <- read.csv("Meta_train.csv")
data2 <- read.csv("Meta_test.csv")

# Merge datasets by Traits column
df <- merge(data1, data2, by = "Traits")

# Load plotting libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# Preprocess data for 9 different volcano plots
for (i in 1:9) {
  # Define column names for each iteration
  HR_col <- paste0("HR_", i, ".x")          # Hazard Ratio from training data
  HR_col1 <- paste0("HR_", i, ".y")         # Hazard Ratio from test data
  FDR_col <- paste0("FDR_Pvalue_", i, ".x") # FDR-adjusted p-value from training
  P_col <- paste0("Pvalue_", i, ".y")       # Raw p-value from test data
  LogP_col <- paste0("LogFDR", i)           # Will store -log10(FDR)
  color_group_col <- paste0("color_group", i) # Will store color categories
  
  # Calculate -log10(FDR) for volcano plot y-axis
  df[[LogP_col]] <- -log10(df[[FDR_col]])
  
  # Handle extreme values (0 or >100 becomes 90 for visualization)
  df[[LogP_col]][df[[LogP_col]] == 0 | df[[LogP_col]] > 100] <- 90
  
  # Create color categories based on significance and effect direction:
  # - Deep Red: Significant in both training and test (HR > 1)
  # - Light Red: Significant in training only (HR > 1)
  # - Deep Blue: Significant in both training and test (HR < 1)
  # - Light Blue: Significant in training only (HR < 1)
  # - Gray: Not significant
  df[[color_group_col]] <- case_when(
    df[[HR_col]] > 1 & df[[HR_col1]] > 1 & df[[FDR_col]] < 0.05 & df[[P_col]] < 0.05 ~ "Deep Red",
    df[[HR_col]] > 1 & df[[HR_col1]] > 1 & df[[FDR_col]] < 0.05 & df[[P_col]] > 0.05 ~ "Light Red",
    df[[FDR_col]] >= 0.05 ~ "Gray",
    df[[HR_col]] < 1 & df[[HR_col1]] < 1 & df[[FDR_col]] < 0.05 & df[[P_col]] < 0.05 ~ "Deep Blue",
    df[[HR_col]] < 1 & df[[HR_col1]] < 1 & df[[FDR_col]] < 0.05 & df[[P_col]] > 0.05 ~ "Light Blue",
    TRUE ~ "Gray"
  )
}

# Function to create a single volcano plot
plot_volcano <- function(df, HR_col, LogP_col, color_group_col) {
  
  # Identify top significant points for labeling:
  # - Top 2 positive associations (HR > 1)
  significant_points_pos <- df %>%
    filter(.data[[HR_col]] > 1) %>%
    arrange(desc(.data[[HR_col]]), desc(.data[[LogP_col]])) %>%
    head(2)
  
  # - Top 2 negative associations (HR < 1)
  significant_points_neg <- df %>%
    filter(.data[[HR_col]] < 1) %>%
    arrange(.data[[HR_col]], desc(.data[[LogP_col]])) %>%
    head(2)
  
  # Combine significant points
  significant_points <- bind_rows(significant_points_pos, significant_points_neg)
  
  # Create the volcano plot
  ggplot(df, aes(x = .data[[HR_col]], y = .data[[LogP_col]], color = .data[[color_group_col]])) +
    ylim(0, 100) +                              # Set y-axis limits
    xlim(0.5, 2) +                              # Set x-axis limits
    geom_point(alpha = 0.6,                     # Semi-transparent points
               size = ifelse(df[[color_group_col]] %in% c("Deep Red", "Deep Blue"), 3, 2)) +  # Size by significance
    scale_color_manual(                         # Custom color scheme
      values = c("Deep Red" = "#8C0000", 
                 "Light Red" = "#C03636",
                 "Deep Blue" = "#00447C", 
                 "Light Blue" = "#4B7299", 
                 "Gray" = "gray")) +
    geom_hline(yintercept = -log10(0.05),       # Significance threshold line
               linetype = "dashed", 
               color = "gray") +
    labs(x = "Hazard Ratio (HR)",               # Axis labels
         y = "-Log10(FDR-adjusted p-value)") +
    geom_text(data = significant_points,        # Label top points
              aes(label = Traits), 
              vjust = -1, 
              color = "black", 
              position = position_jitter(width = 0.1, height = 0.1)) +  
    theme_minimal() +                           # Clean theme
    theme(panel.border = element_rect(colour = "gray", fill = NA, size = 0.5),
          legend.position = "none")             # Remove legend
}

# Generate all 9 volcano plots
plots <- lapply(1:9, function(i) {
  HR_col <- paste0("HR_", i, ".x")
  LogP_col <- paste0("LogFDR", i)
  color_group_col <- paste0("color_group", i)
  plot_volcano(df, HR_col, LogP_col, color_group_col)
})

# Arrange all 9 plots in a 3x3 grid
library(gridExtra)
grid.arrange(grobs = plots, nrow = 3)