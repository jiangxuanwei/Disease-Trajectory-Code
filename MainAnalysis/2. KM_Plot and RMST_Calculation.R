# ============================================================
# KM Plot and RMST Calculation
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

# Load required libraries
library(dplyr)       # For data manipulation
library(survival)    # For survival analysis
library(survminer)   # For survival visualization
library(survRM2)     # For restricted mean survival time analysis

# Load and merge datasets ------------------------------------------------------

# Read primary dataset with disease transition information
data1 <- read.csv("CMD_CA.csv")

# Read covariates data (previously imputed)
covariate <- read.csv("covariates_impute.csv")

track <- read.csv("track.csv")

# Merge all datasets by participant ID
data <- merge(data1, covariate, by = "n_eid")
data <- merge(data, track, by = "n_eid")

# Create composite cardiovascular disease (CD) event ---------------------------
data$event_CD <- 0
data$event_CD[data$event_dia == 1 | data$event_CVD == 1] <- 1  # Diabetes or CVD event
data$time_CD <- pmin(data$time_dia, data$time_CVD, na.rm = TRUE)  # Time to first CD event

# Define disease transition states and survival times --------------------------

# Status coding:
# 0 = Healthy (no transitions)
# 1 = Only CMD (transition 1 only)
# 2 = Only CA (transition 2 only)
# 3 = CMD → CMD+CA (transitions 1 & 4)
# 4 = CA → CA+CMD (transitions 2 & 6)

# Healthy participants (no transitions)
data$status[data$trans1 == 0 & data$trans2 == 0] <- 0
data$survtime[data$trans1 == 0 & data$trans2 == 0] <- data$death2_time[data$trans1 == 0 & data$trans2 == 0]

# Only CMD transition
data$status[data$trans1 == 1 & data$trans2 == 0] <- 1
data$survtime[data$trans1 == 1 & data$trans2 == 0] <- data$death2_time[data$trans1 == 1 & data$trans2 == 0] - data$time_CD[data$trans1 == 1 & data$trans2 == 0]

# Only CA transition
data$status[data$trans1 == 0 & data$trans2 == 1] <- 2
data$survtime[data$trans1 == 0 & data$trans2 == 1] <- data$death2_time[data$trans1 == 0 & data$trans2 == 1] - data$time_CA[data$trans1 == 0 & data$trans2 == 1]

# CMD → CMD+CA transition
data$status[data$trans1 == 1 & data$trans4 == 1] <- 3
data$survtime[data$trans1 == 1 & data$trans4 == 1] <- data$death2_time[data$trans1 == 1 & data$trans4 == 1] - data$time_CA[data$trans1 == 1 & data$trans4 == 1]

# CA → CA+CMD transition
data$status[data$trans2 == 1 & data$trans6 == 1] <- 4
data$survtime[data$trans2 == 1 & data$trans6 == 1] <- data$death2_time[data$trans2 == 1 & data$trans6 == 1] - data$time_CD[data$trans2 == 1 & data$trans6 == 1]

# Convert status to factor for proper grouping
data$status <- as.factor(data$status)

# Create cumulative incidence plot ---------------------------------------------

# Fit survival curves by transition status
fit0 <- survfit(Surv(survtime, death2_icd10) ~ status, data = data)

# Generate cumulative incidence plot
p <- ggsurvplot(fit0, data = data, size = 1,
                fun = "event",  # Show cumulative events rather than survival
                palette = c("grey", "#A8585B", "#DB946D", "#7C739E", "#1A4789"),  # Custom colors
                risk.table = FALSE,  # Don't show risk table
                pval = " ",  # Empty space where p-value would go
                legend.title = element_text(''),  # No legend title
                conf.int = TRUE,  # Show confidence intervals
                censor = FALSE,  # Don't show censoring marks
                legend.labs = c("Health", "Only CMD", "Only CA", "CMD-CA", "CA-CMD"),  # Group labels
                xlim = c(0, 15),  # X-axis limit (15 years)
                ylim = c(0, 0.6),  # Y-axis limit (60% cumulative incidence)
                break.time.by = 5,  # X-axis ticks every 5 years
                ylab = "Cumulative incidence of all-cause mortality",
                xlab = "Time since diagnosis (years)",
                risk.table.y.text = FALSE) 

# Display the plot
p

# Restricted Mean Survival Time (RMST) Analysis --------------------------------

# Prepare data for CA+CMD vs Healthy comparison,comparison for other states also follow this pipelines
data1 <- data %>% 
  filter(status %in% c(0, 4) & !is.na(survtime)) %>%  # Keep only healthy and CA+CMD
  mutate(status = ifelse(status == 4, 1, 0))  # Recode status as 0/1 for analysis

# Calculate RMST difference at 15 years
rmst_result <- rmst2(time = data1$survtime, 
                     status = data1$death2_icd10,
                     arm = data1$status, 
                     tau = 15)  # Time horizon of 15 years
