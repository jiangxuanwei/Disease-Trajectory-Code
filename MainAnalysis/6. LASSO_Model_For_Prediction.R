# ============================================================
# Prediction model construction
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

# Load required packages
library(glmnet)       # For LASSO regression
library(dplyr)        # For data manipulation
library(survival)     # For survival analysis

# Load and prepare protein data

# Read training and test protein results
datatrans <- read.csv("Prot_train.csv")
data2 <- read.csv("Prot_test.csv")

# Merge training and test results by protein traits
data <- merge(datatrans, data2, by="Traits")

# Filter significant proteins for each transition (FDR < 0.05 in both training and test)
for (i in 1:9) {
  FDR_col <- paste0("FDR_Pvalue_", i, ".x")  # Training FDR column
  Pvalue_col <- paste0("Pvalue_", i, ".y")   # Test p-value column
  
  # Get proteins significant in both datasets
  Prot_filtered <- data %>%
    filter(get(FDR_col) < 0.05 & get(Pvalue_col) < 0.05) %>%
    pull(Traits)
  
  # Store filtered protein names for each transition
  assign(paste0("id", i), Prot_filtered)
}

# Prepare main dataset
datatrans <- read.csv("CMD_CA.csv")

# Read and standardize protein measurements
dataProt <- read.csv("protein_UKB_filled.csv")
dataProt[, 2:2921] <- scale(dataProt[, 2:2921])  # Standardize all proteins

# Read covariates
corv <- read.csv("covariates_impute.csv")

# Merge all data
datatrans <- merge(datatrans, dataProt, by="n_eid")
datatrans <- merge(datatrans, corv, by="n_eid")

# Create outcome variables for each transition

# Transition 1: Healthy → CMD
datatrans$outcome1 <- 0
datatrans$outcome1[datatrans$trans1==1 & datatrans$trans4==0 & datatrans$trans5==0] <- 1

# Transition 2: Healthy → Cancer
datatrans$outcome2 <- 0
datatrans$outcome2[datatrans$trans2==1 & datatrans$trans6==0 & datatrans$trans7==0] <- 1

# Transition 3: Healthy → Death
datatrans$outcome3 <- 0
datatrans$outcome3[datatrans$trans3==1] <- 1

# Transition 4: CMD → CMD+Cancer
datatrans$outcome4 <- 0
datatrans$outcome4[datatrans$trans4==1 & datatrans$trans8==0] <- 1

# Transition 5: CMD → Death
datatrans$outcome5 <- 0
datatrans$outcome5[datatrans$trans5==1] <- 1

# Transition 6: Cancer → CMD+Cancer
datatrans$outcome6 <- 0
datatrans$outcome6[datatrans$trans6==1 & datatrans$trans9==0] <- 1

# Transition 7: Cancer → Death
datatrans$outcome7 <- 0
datatrans$outcome7[datatrans$trans7==1] <- 1

# Transition 8: CMD+Cancer → Death
datatrans$outcome8 <- 0
datatrans$outcome8[datatrans$trans8==1] <- 1

# Transition 9: CMD → Death
datatrans$outcome9 <- 0
datatrans$outcome9[datatrans$trans9==1] <- 1

# Split data into training and test sets
train.data <- datatrans %>% filter(birthplace == 1)        # Training set (birthplace=1)
test.data <- datatrans %>% filter(is.na(birthplace) | birthplace != 1)  # Test set

### PART 1: LASSO Cox regression for each transition -------------------------------------
for (i in 1:9) {
  # Create predictor matrix (Age, Sex + significant proteins)
  text <- paste0("x <- as.matrix(train.data %>% select(all_of(c('Age','Sex', id", i, "))))")
  eval(parse(text = text))
  
  # Create survival outcome
  y <- eval(parse(text = paste0("Surv(train.data$time", i, "to, train.data$outcome", i, ")")))
  
  # Cross-validate to find optimal lambda
  cv <- cv.glmnet(x, y, family = 'cox', alpha = 1)  # alpha=1 for LASSO
  lambda_best <- cv$lambda.min  # Lambda with minimum cross-validation error
  
  # Fit final LASSO model with optimal lambda
  model <- glmnet(x, y, family = 'cox', alpha = 1, lambda = lambda_best)
  
  # Extract and store coefficients
  coefs <- coef(model)[, ]
  assign(paste0("coef", i), coefs)
}

# Save LASSO coefficients for each transition
write.csv(coef1, "coef_Prot1.csv")
write.csv(coef2, "coef_Prot2.csv") 
write.csv(coef3, "coef_Prot3.csv")
write.csv(coef4, "coef_Prot4.csv")
write.csv(coef5, "coef_Prot5.csv")
write.csv(coef6, "coef_Prot6.csv")
write.csv(coef7, "coef_Prot7.csv")
write.csv(coef8, "coef_Prot8.csv")
write.csv(coef9, "coef_Prot9.csv")

### Genomics and metabolomics data follow the same analysis pipeline as described above

### PART 2: Performance in the testing dataset -------------------------------------
# Function to format concordance index (C-index) with 95% CI
FUNC = function(model){
  c = round(summary(model)$concordance[1], 2)  # C-index
  c.lwr = round(summary(model)$concordance[1] - 1.96*summary(model)$concordance[2], 2)  # Lower CI
  c.upr = round(summary(model)$concordance[1] + 1.96*summary(model)$concordance[2], 2)  # Upper CI
  paste0(c, " (", c.lwr, ", ", c.upr, ")")    # Formatted string
}

# Function to format change in concordance index between models
FUNChange = function(model0, model1){
  modelcompare <- CsChange(model0, model1, data=test.data, nb=100)  # Compare models
  change = round(modelcompare[[1]]$change, 2)      # Change in C-index
  change.lwr = round(modelcompare[[1]]$low, 2)     # Lower CI
  change.upr = round(modelcompare[[1]]$up, 2)      # Upper CI
  paste0(change, " (", change.lwr, ", ", change.upr, ")")  # Formatted string
}

# Initialize data frames for storing results
cstatisticomics <- data.frame(matrix(ncol=9, nrow=1))  # For protein model C-indices
cstatisticbase <- data.frame(matrix(ncol=9, nrow=1))   # For baseline model C-indices
cchange <- data.frame(matrix(ncol=9, nrow=1))          # For C-index changes

# Calculate protein scores and evaluate models ---------------------------------
for (i in 1:9) {
  # Create protein score for current transition
  ids <- paste0("id", i)
  text <- paste0("Metaset <- datatrans %>% select(all_of(c('Age','Sex', id", i, ")))")
  eval(parse(text=text))
  Metaset_matrix <- as.matrix(Metaset)
  
  # Get coefficients for current transition
  text1 <- paste0("current_coefs <- coef", i)
  eval(parse(text=text1))
  
  # Calculate protein score (linear combination of selected proteins)
  datatrans[[paste0("Metascore", i)]] <- Metaset_matrix %*% current_coefs[,2]
  
  # Create test set (non-birthplace=1 samples)
  test.data <- datatrans %>% filter(is.na(birthplace) | birthplace != 1)
  
  # Fit baseline model (Age + Sex only)
  text1 <- paste0("model0 <- coxph(Surv(time", i, "to, outcome", i, ") ~ Age + Sex, data = test.data)")
  
  # Fit protein model (protein score only)
  text2 <- paste0("model1 <- coxph(Surv(time", i, "to, outcome", i, ") ~ Metascore", i, ", data = test.data)")
  
  eval(parse(text=text1))
  eval(parse(text=text2))
  
  # Store results
  cstatisticbase[, i] <- FUNC(model0)    # Baseline model performance
  cstatisticomics[, i] <- FUNC(model1)   # Protein model performance
  cchange[, i] <- FUNChange(model0, model1)  # Improvement from adding proteins
}

# Combine and save results 
result <- rbind(cstatisticbase, cstatisticomics, cchange)

# Prepare and save output dataset with protein scores
dataout <- datatrans[, c("n_eid", 
                         paste0("trans", 1:9),
                         paste0("time", 1:9, "from"),
                         paste0("time", 1:9, "to"),
                         paste0("outcome", 1:9),
                         paste0("Metascore", 1:9))]

write.csv(dataout, "ProtScore.csv")
write.csv(result, "Cindex Prot.csv")

### Genomics and metabolomics data follow the same analysis pipeline as described above