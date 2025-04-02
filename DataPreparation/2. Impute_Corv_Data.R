# ============================================================
# Imputation
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

# Perform multiple imputation with 5 datasets ----------------------------------
data <- read.csv("covariates.csv")
library(dplyr)

# Convert categorical variables to factors -------------------------------------
data$Sex <- as.factor(data$Sex)
data$Ethnicity <- as.factor(data$Ethnicity)
data$Edu <- as.factor(data$Edu)
data$Employed <- as.factor(data$Employed)
data$Smoke <- as.factor(data$Smoke)
data$Drink <- as.factor(data$Drink)
data$T2dhis <- as.factor(data$T2dhis)
data$Cvdhis <- as.factor(data$Cvdhis)
data$hiscancer <- as.factor(data$hiscancer)
data$centre <- as.factor(data$centre)
data$birthplace <- as.factor(data$birthplace)
data$tumorstage <- as.factor(data$tumorstage)
data$cancer_screen <- as.factor(data$cancer_screen)
data$immune_treat <- as.factor(data$immune_treat)
data$cell_treat <- as.factor(data$cell_treat)
data$CMD_treat <- as.factor(data$CMD_treat)
library(mice)
# Select variables to impute (same as before)
data_impute <- data[, c("n_eid", "Ethnicity", "Edu", "Employed", "Tdi","Smoke", "Drink", "Mets",
                        "Sleeptime", "cancer_screen", "glucose", "ldl", "sbp", "dbp", "cho", "hdl")]

# Run multiple imputation (5 imputations)
imp <- mice(data_impute, m = 5, maxit = 5, defaultMethod = c("norm", 'logreg', 'polyreg', 'polr'))

# Check convergence
plot(imp)

# Create a list to store the 5 completed datasets
completed_data <- list()
for (i in 1:5) {
  completed_data[[i]] <- complete(imp, i)
}

### Create a pooled dataset (simple average for continuous variables)
pooled_data <- completed_data[[1]]  # Start with first imputation
numeric_vars <- sapply(pooled_data, is.numeric)

# Average across imputations for numeric variables
for (var in names(pooled_data)[numeric_vars]) {
  temp <- sapply(completed_data, function(x) x[[var]])
  pooled_data[[var]] <- rowMeans(temp)
}

# For factors, take the mode (most common value)
factor_vars <- sapply(pooled_data, is.factor)
for (var in names(pooled_data)[factor_vars]) {
  temp <- sapply(completed_data, function(x) as.character(x[[var]]))
  pooled_data[[var]] <- apply(temp, 1, function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  })
  pooled_data[[var]] <- as.factor(pooled_data[[var]])
}

# Merge with non-imputed variables (same as before)
data1 <- data[,c("n_eid","Age","Sex","Dietscore","BMI","T2dhis","Cvdhis","hiscancer","centre","birthplace",
                 "tumorstage","cancer_treat","CMD_treat")]
data_imp_pooled <- merge(data1, pooled_data, by="n_eid")
