# ============================================================
# Principal Component
# Author: JXW, YGR
# Last Modified: Feb 2025
# ============================================================

### 1. At the individual level--------------------------------
### Proteomics-----
# Load the required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Read in the necessary data files
data <- read.csv("CMD_CA.csv")  # Main data containing transitions and death time
protein <- read.csv("protein_UKB_filled.csv")  # Proteomics data

# Merge the data on "n_eid" from the CMD_CA and protein datasets
data1 <- merge(data, protein, by.x = "n_eid", by.y = "eid")

# Read the significant protein markers list
signpro <- read.csv("Protsign.csv")
signpros <- signpro$meta  # Extract protein markers of interest

# Prepare the PCA dataset by assigning status and calculating survival time based on disease transitions
PCAset <- data1

# Assign disease status and calculate survival time for each transition state
PCAset$status[PCAset$trans1 == 0 & PCAset$trans2 == 0] <- 0
PCAset$survtime[PCAset$trans1 == 0 & PCAset$trans2 == 0] <- PCAset$death2_time[PCAset$trans1 == 0 & PCAset$trans2 == 0]

PCAset$status[PCAset$trans1 == 1 & PCAset$trans2 == 0] <- 1
PCAset$survtime[PCAset$trans1 == 1 & PCAset$trans2 == 0] <- PCAset$death2_time[PCAset$trans1 == 1 & PCAset$trans2 == 0] - PCAset$time_CD[PCAset$trans1 == 1 & PCAset$trans2 == 0]

PCAset$status[PCAset$trans1 == 0 & PCAset$trans2 == 1] <- 2
PCAset$survtime[PCAset$trans1 == 0 & PCAset$trans2 == 1] <- PCAset$death2_time[PCAset$trans1 == 0 & PCAset$trans2 == 1] - PCAset$time_CA[PCAset$trans1 == 0 & PCAset$trans2 == 1]

PCAset$status[PCAset$trans1 == 1 & PCAset$trans4 == 1] <- 3
PCAset$survtime[PCAset$trans1 == 1 & PCAset$trans4 == 1] <- PCAset$death2_time[PCAset$trans1 == 1 & PCAset$trans4 == 1] - PCAset$time_CA[PCAset$trans1 == 1 & PCAset$trans4 == 1]

PCAset$status[PCAset$trans2 == 1 & PCAset$trans6 == 1] <- 4
PCAset$survtime[PCAset$trans2 == 1 & PCAset$trans6 == 1] <- PCAset$death2_time[PCAset$trans2 == 1 & PCAset$trans6 == 1] - PCAset$time_CD[PCAset$trans2 == 1 & PCAset$trans6 == 1]

# Convert status to a factor type
PCAset$status <- as.factor(PCAset$status)

# Extract only the significant proteins that exist in the PCAset data
existing_proteins <- intersect(signpros, colnames(PCAset))

# Prepare the PCA dataset with status and relevant proteins
pca_data <- PCAset %>%
  select(status, all_of(existing_proteins))

# Convert status to factor with labels for better readability
pca_data$status <- factor(pca_data$status, 
                          levels = 0:4,
                          labels = paste("Status", 0:4))

# Perform PCA on the selected proteins, centering and scaling the data
pca_result <- prcomp(pca_data[, existing_proteins], 
                     center = TRUE, 
                     scale. = TRUE)

# Extract PCA scores (PC1, PC2, etc.)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$status <- pca_data$status

# Calculate the percentage of variance explained by each principal component
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create the PCA plot
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = status)) +
  geom_point(size = 3, alpha = 0.7) +  # Plot the points with transparency
  stat_ellipse(level = 0.95) +  # Add 95% confidence ellipses around points
  scale_color_manual(values = c("Status 0" = "grey", 
                                "Status 1" = "#A8585B", 
                                "Status 2" = "#DB946D", 
                                "Status 3" = "#7C739E", 
                                "Status 4" = "#1A4789")) +  # Custom colors for each status
  labs(
    title = "PCA of Significant Proteins by Disease Status",  # Title for the plot
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),  # Label for x-axis (PC1)
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),  # Label for y-axis (PC2)
    color = "Disease Status"  # Legend label for the color mapping
  ) +
  theme_bw() +  # Use a black-and-white theme for the plot
  theme(
    legend.position = "right",  # Position the legend on the right
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center the title and make it bold
    axis.title = element_text(face = "bold")  # Make axis titles bold
  )


### 2. At the age group level--------------------------------
### Proteomics-----
# Read the proteomics data (with imputed missing values)
datameta <- read.csv("protein_UKB_filled.csv")

# Standardize the protein expression matrix (from column 2 to column 2921) by Z-score
datameta[, 2:2921] <- scale(datameta[, 2:2921])

# Read the transition data (which tracks disease progression and event status)
trans <- read.csv("CMD_cA.csv")

# Merge the proteomics data with transition data using "n_eid" as the key
data <- merge(datameta, trans, by = "n_eid")

# Read covariate data
corv <- read.csv("covariates_impute.csv")

# Merge the covariate data with the previously merged dataset
data <- merge(data, corv, by = "n_eid")

# Read the significant protein markers list
signprot <- read.csv("proteinsign.csv")

# Create the status variable based on transition events
data$status[data$trans1 == 0 & data$trans2 == 0] <- 0
data$status[data$trans1 == 1 & data$trans2 == 0] <- 1
data$status[data$trans1 == 0 & data$trans2 == 1] <- 2
data$status[data$trans1 == 1 & data$trans4 == 1] <- 3
data$status[data$trans2 == 1 & data$trans6 == 1] <- 4

# Categorize individuals into age groups
data$age_group <- cut(data$Age, 
                      breaks = c(-Inf, 50, 55, 60, 65, Inf),  
                      labels = c("<50", "50-55", "55-60", "60-65", ">65"))

# Combine the status and age group variables to create a new interaction variable
data$status_age_group <- interaction(data$status, data$age_group)

# Select proteins from the list of significant proteins that are available in the dataset
protein_vars <- signprot$proteins
available_proteins <- protein_vars[protein_vars %in% colnames(data)]

# Aggregate the data to calculate the mean of each protein for each status and age group combination
library(dplyr)
protein_means <- aggregate(data[, available_proteins], 
                           by = list(status_age_group = data$status_age_group), 
                           FUN = function(x) mean(x, na.rm = TRUE))

# Assign the status labels for each group
protein_means$status <- rep(c(1:5), 5)

# Remove the first column (status_age_group) and the 91st column (which is likely a non-relevant column) from the aggregated data
protein_means_data <- protein_means[, c(-1, -91)]

# Standardize the aggregated protein data
protein_means_data_scaled <- scale(protein_means_data)

# Perform Principal Component Analysis (PCA) on the scaled data
pca_result <- prcomp(protein_means_data_scaled, center = TRUE, scale. = TRUE)

# Extract the PCA scores
pca_scores <- data.frame(pca_result$x)
pca_scores$status_age_group <- protein_means$status_age_group
pca_scores$status <- protein_means$status

# Define the colors to represent different statuses in the PCA plot
status_colors <- c("grey", "#A8585B", "#DB946D", "#7C739E", "#1A4789")

# Create the PCA plot using ggplot2
ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(status))) +
  geom_point(alpha = 0.7, size = 3) +  # Plot the points with transparency
  scale_color_manual(values = status_colors) +  # Customize the colors for different statuses
  stat_ellipse(aes(x = PC1, y = PC2, color = as.factor(status)), 
               level = 0.95,  # Add ellipses around points based on status with a confidence level of 95%
               linewidth = 0.5, 
               alpha = 1) + 
  geom_text(aes(label = status_age_group), vjust = -1, size = 3) +  # Add labels to the points with the status_age_group
  theme_minimal() +  # Use a minimal theme
  labs(
    title = "PCA of Protein Means by Status and Age Group",  # Plot title
    x = "Principal Component 1",  # X-axis label
    y = "Principal Component 2",  # Y-axis label
    color = "Status & Age Group"  # Legend title for the color
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Customize plot title style
    legend.title = element_text(size = 12),  # Customize legend title style
    legend.text = element_text(size = 10)  # Customize legend text style
  )

### Metabolomics data also follow the above pipelines