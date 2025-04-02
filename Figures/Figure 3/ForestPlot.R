# ============================================================
# Forest Plot
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

# Load necessary libraries
library(data.table)   # For reading and handling large data files
library(forestploter) # For creating forest plots
library(grid)         # For graphical parameters
library(dplyr)        # For data manipulation
library(ggplot2)      # For general data visualization
library(ggpubr)       # For enhanced ggplot2 functionalities
library(ggsci)        # For additional color palettes

# Load the data from the CSV file using fread (efficient for large datasets)
data <- fread("plot.csv", header = TRUE)

# Create an empty column with spaces for formatting purposes in the plot
data$` ` <- paste(rep(" ", 20), collapse = " ")  

# Generate a formatted string for the C-index with its 95% confidence intervals
data$Cindex111 <- ifelse(is.na(data$Cindex), "", 
                         sprintf('%.2f (%.2f - %.2f)', data$Cindex, data$LCI, data$UCI))

# Create a custom theme for the forest plot
tm <- forest_theme(base_size = 10,               # Set base font size
                   ci_pch = 20,                  # Circle symbol for confidence interval
                   ci_col = "black",             # Confidence interval color
                   ci_lty = 1,                   # Confidence interval line type
                   ci_lwd = 1,                   # Confidence interval line width
                   ci_Theight = 0,               # Height of CI labels
                   refline_gp = gpar(lwd = 1.5, lty = "dashed", col = "black"), # Reference line styling
                   summary_fill = "black",       # Fill color for summary box
                   summary_col = "black",        # Color of summary text
                   footnote_gp = gpar(cex = 0.7, fontface = "italic", col = "blue"))  # Footnote styling

# Create the forest plot using the forest function from the forestploter package
p <- forest(data[,c(1,6,5)],                      # Select relevant columns from data (name, C-index, LCI, UCI)
            est = data$Cindex,                    # Effect size (C-index)
            lower = data$LCI,                     # Lower bound of the confidence interval
            upper = data$UCI,                     # Upper bound of the confidence interval
            sizes = 0.8,                          # Adjust the size of the plot markers
            ci_column = 3,                        # Column containing the confidence interval
            xlim = c(0.5, 1),                     # Set the x-axis limits (C-index range)
            ticks_at = c(0.5, 0.6, 0.7, 0.8, 0.9, 1), # Set ticks on x-axis
            theme = tm)                           # Apply the custom theme

# Print the plot to the R console
p

# Save the plot as a PDF file with specified dimensions (width = 5.5, height = 5.5)
ggsave(paste("forest.pdf", sep = ""), p, width = 5.5, height = 5.5)
