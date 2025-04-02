# ============================================================
# Bubble Plot
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

# Load required libraries
library(ggplot2)    # For data visualization
library(dplyr)      # For data manipulation
library(tidyr)      # For data reshaping
library(tidyverse)  # Collection of packages for data science

### 1. Genomics (PRS) Bubble Plot ###
# Load polygenic risk score (PRS) coefficients data
df <- read.csv("PRS_coef_merge.csv", 
               na.strings = c("", " "),  # Treat empty strings as NA
               stringsAsFactors = FALSE)

# Reshape data from wide to long format for plotting
df <- df %>%
  pivot_longer(cols = everything(),  # Transform all columns
               names_to = "Variable",  # Coefficient names will go here
               values_to = "ID") %>%  # Protein IDs will go here
  drop_na()  # Remove rows with missing values

# Count how many times each protein appears across coefficients
prot_counts <- df %>%
  group_by(ID) %>%
  summarise(Count = n())  # Count occurrences per protein

# Merge counts back with original data
df <- merge(df, prot_counts, by = "ID")

# Filter out any remaining NA values
df <- df %>% filter(!is.na(ID))

# Order factors by count (most frequent proteins first)
df$ID <- factor(df$ID, levels = prot_counts$ID[order(prot_counts$Count, decreasing = TRUE)])
df$Variable <- factor(df$Variable, levels = paste0("coef", 1:9))  # Force coefficient order

# Create bubble plot for genomic markers
ggplot(df, aes(x = ID, y = Variable, color = Count)) +
  geom_point(size = 8, alpha = 0.7) +  # Size and transparency for bubbles
  theme_minimal() +  # Clean background
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(title = "Genomic Markers Association",
       x = "Protein",
       y = "Transition States") +
  scale_y_discrete(drop = FALSE) +  # Keep all coefficient levels even if empty
  scale_color_gradient(high = "#064F77", low = "#67889B")  # Blue color gradient

### 2. Metabolomics Bubble Plot ###
# Load metabolomics coefficients data
df <- read.csv("Meta_coef_merge.csv", 
               na.strings = c("", " "), 
               stringsAsFactors = FALSE)

# Load metabolite name mappings
marker_name <- read.csv("Meta markers.csv")
marker_suo <- read.csv("Meta_abbre.csv")

# Data preprocessing pipeline:
df_1 <- df %>%
  pivot_longer(cols = everything(), 
               names_to = "Variable", 
               values_to = "ID") %>%
  drop_na()

# Clean metabolite IDs by removing prefixes/suffixes
df_2 <- df_1 %>%
  mutate(ID = str_replace(ID, "n_", "")) %>%  # Remove 'n_' prefix
  mutate(ID = str_replace(ID, "_0_0", ""))    # Remove '_0_0' suffix

# Merge with name mappings
df_3 <- merge(df_2, marker_name, by = "ID")
df_4 <- merge(df_3, marker_suo, by = "Name")

# Count metabolites and filter for those appearing ≥3 times
metab_counts <- df_4 %>%
  group_by(Metasuo) %>%  # Group by metabolite abbreviation
  summarise(Count = n()) %>%
  filter(Count >= 3)     # Minimum occurrence threshold

# Merge counts back and prepare factors
df_4 <- merge(df_4, metab_counts, by = "Metasuo")
df_4$Metasuo <- factor(df_4$Metasuo, levels = metab_counts$Metasuo[order(metab_counts$Count, decreasing = TRUE)])
df_4$Variable <- factor(df_4$Variable, levels = paste0("coef", 1:9))

# Create metabolomics bubble plot
ggplot(df_4, aes(x = Metasuo, y = Variable, color = Count)) +
  geom_point(size = 8, alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Metabolomic Markers Association",
       x = "Metabolite",
       y = "Transition States") +
  scale_y_discrete(drop = FALSE) +
  scale_color_gradient(high = "#FF7F00", low = "#EFB454")  # Orange color gradient

### 3. Proteomics Bubble Plot ###
# Load proteomics coefficients data
df <- read.csv("Prot_coef_merge.csv", 
               na.strings = c("", " "), 
               stringsAsFactors = FALSE)

# Process similarly to genomic data
df <- df %>%
  pivot_longer(cols = everything(), 
               names_to = "Variable", 
               values_to = "ID") %>%
  drop_na()

# Count and filter proteins appearing >3 times
prot_counts <- df %>%
  group_by(ID) %>%
  summarise(Count = n()) %>%
  filter(Count > 3)  # More stringent threshold for proteomics

df <- merge(df, prot_counts, by = "ID")
df <- df %>% filter(!is.na(ID))

# Order factors
df$ID <- factor(df$ID, levels = prot_counts$ID[order(prot_counts$Count, decreasing = TRUE)])
df$Variable <- factor(df$Variable, levels = paste0("coef", 1:9))

# Create proteomics bubble plot
ggplot(df, aes(x = ID, y = Variable, color = Count)) +
  geom_point(size = 8, alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Proteomic Markers Association",
       x = "Protein",
       y = "Transition States") +
  scale_y_discrete(drop = FALSE) +
  scale_color_gradient(high = "#B33D29", low = "#EA9C92")  # Red color gradient