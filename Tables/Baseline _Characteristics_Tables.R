# ============================================================
# Baseline Characteristics
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

# Load necessary datasets for the analysis
data1 <- read.csv("CMD_CA.csv")
geno <- read.csv("geno_dataset.csv")
meta <- read.csv("meta_dataset.csv")
prot <- read.csv("prote_dataset.csv")
data2 <- read.csv("covariates_impute.csv")

# Merge the datasets by common identifier 'n_eid'
data <- merge(data1, data2, by = "n_eid")
data <- merge(geno, data, by = "n_eid")
data <- merge(meta, data, by = "n_eid")
data <- merge(prot, data, by = "n_eid")

# Convert specific variables to factor type
factor_vars <- c("Sex", "Ethnicity", "Edu", "Employed", 
                 "Smoke", "Drink", "T2dhis", "Cvdhis", 
                 "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen")

data <- data %>%
  mutate(across(all_of(factor_vars), as.factor))  # Convert specified variables to factor type

# Create a new variable `status` based on the values of transition variables
data <- data %>%
  mutate(status = case_when(
    trans1 == 0 & trans2 == 0 ~ 0,  # No transitions (status 0)
    trans1 == 1 & trans4 == 0 ~ 1,  # Transition 1 occurred (status 1)
    trans2 == 1 & trans6 == 0 ~ 2,  # Transition 2 occurred (status 2)
    trans4 == 1 | trans6 == 1 ~ 3,  # Transition 4 or 6 occurred (status 3)
    TRUE ~ 4  # Default case (status 4)
  ))

# Generate a descriptive summary table based on the status variable
library(table1)
table1(~ Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore + Sleeptime + 
         BMI + T2dhis + Cvdhis + hiscancer + cancer_treat + CMD_treat + cancer_screen | status, data)

### 2. HRS

# Load data for HRS analysis
data1 <- read.csv("data_no_impute.csv")

# Handle BMI outliers by setting BMI values outside of a certain range to NA
data1$BMI[data1$BMI < 5 | data1$BMI > 50] <- NA

# Reload imputed data and merge with other relevant data
data1 <- read.csv("data_impute.csv")
data2 <- read.csv("CMD_CA_trans_HRS.csv")
data <- merge(data1, data2, by = "HHIDPN")

# Convert specified variables to factor type
factor_vars <- c("Sex", "Edu", "Employed", "Income", "Smoke", "Drink", "Mets", "CMD_treat")
data <- data %>%
  mutate(across(all_of(factor_vars), as.factor))

# Create a new variable `status` based on transition variables
data <- data %>%
  mutate(status = case_when(
    trans1 == 0 & trans2 == 0 ~ 0,  # No transitions (status 0)
    trans1 == 1 & trans4 == 0 ~ 1,  # Transition 1 occurred (status 1)
    trans2 == 1 & trans6 == 0 ~ 2,  # Transition 2 occurred (status 2)
    trans4 == 1 | trans6 == 1 ~ 3,  # Transition 4 or 6 occurred (status 3)
    TRUE ~ 4  # Default case (status 4)
  ))

# Generate a descriptive summary table based on the status variable
library(table1)
table1(~ Age + Sex + Ethnicity + Edu + Employed + 
         Income + Smoke + Drink + Mets + BMI | status, data = data)

### 3. After PSM

# Load matched survival data
data <- read.csv("match_data.csv")

# Convert specified variables to factor type
factor_vars <- c("Sex", "Ethnicity", "Edu", "Employed", 
                 "Smoke", "Drink", "T2dhis", "Cvdhis", 
                 "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen")

data <- data %>%
  mutate(across(all_of(factor_vars), as.factor))

# Assign status and survival time based on the transition variables
data$status[data$trans1 == 0 & data$trans2 == 0] <- 0  # No transitions (status 0)
data$survtime[data$trans1 == 0 & data$trans2 == 0] <- data$death2_time[data$trans1 == 0 & data$trans2 == 0]

data$status[data$trans1 == 1 & data$trans2 == 0] <- 1  # Transition 1 occurred (status 1)
data$survtime[data$trans1 == 1 & data$trans2 == 0] <- data$death2_time[data$trans1 == 1 & data$trans2 == 0] - data$time_CD[data$trans1 == 1 & data$trans2 == 0]

data$status[data$trans1 == 0 & data$trans2 == 1] <- 2  # Transition 2 occurred (status 2)
data$survtime[data$trans1 == 0 & data$trans2 == 1] <- data$death2_time[data$trans1 == 0 & data$trans2 == 1] - data$time_CA[data$trans1 == 0 & data$trans2 == 1]

data$status[data$trans1 == 1 & data$trans4 == 1] <- 3  # Transition 4 occurred (status 3)
data$survtime[data$trans1 == 1 & data$trans4 == 1] <- data$death2_time[data$trans1 == 1 & data$trans4 == 1] - data$time_CA[data$trans1 == 1 & data$trans4 == 1]

data$status[data$trans2 == 1 & data$trans6 == 1] <- 4  # Transition 6 occurred (status 4)
data$survtime[data$trans2 == 1 & data$trans6 == 1] <- data$death2_time[data$trans2 == 1 & data$trans6 == 1] - data$time_CD[data$trans2 == 1 & data$trans6 == 1]

data$status <- as.factor(data$status)  # Convert status to factor type

# Generate a descriptive summary table based on the status variable
library(table1)
table1(~ Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore + Sleeptime + 
         BMI + T2dhis + Cvdhis + hiscancer + cancer_treat + CMD_treat + cancer_screen | status, data)
