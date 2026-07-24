library(dplyr)
library(survival)
library(glmnet)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

data1 <- read.csv("D:/科研数据/课题数据/1. UKB数据课题/UKB数据常用汇总/Obesity data/BMI.csv")
data2 <- read.csv("D:/科研数据/课题数据/1. UKB数据课题/UKB数据常用汇总/Omics data/protein_UKB_filled.csv")
data3 <- read.csv("D:/科研数据/课题数据/1. UKB数据课题/UKB数据常用汇总/covariates_impute.csv")
track <- read.csv("D:/科研数据/课题数据/1. UKB数据课题/All disease/UBK_disease.csv")

# Log-transform and standardize proteomics data.
#data2[, 2:2921] <- scale(log2(data2[, 2:2921] + 20))

# Merge BMI, proteomics, covariate, and disease follow-up data.
data <- merge(data1, data2, by = "n_eid")
data <- merge(data, data3, by = "n_eid")
data <- merge(data, track, by = "n_eid")

# Fill missing disease survival times with death follow-up time.
data <- data %>%
  mutate(
    across(
      starts_with("survival_time_"),
      ~ ifelse(is.na(.), death2_time, .)
    )
  )

# Define cancer event and follow-up time variables.
cancer_events <- c(
  "event_LungCA", "event_BreastCA", "event_ColorectalCA",
  "event_ProstaticCA", "event_GastroCA", "event_liverCA",
  "event_EsophagealCA", "event_PancreaticCA",
  "event_CervicalCA", "event_ovarianCA"
)

cancer_times <- c(
  "time_LungCA", "time_BreastCA", "time_ColorectalCA",
  "time_ProstaticCA", "time_GastroCA", "time_liverCA",
  "time_EsophagealCA", "time_PancreaticCA",
  "time_CervicalCA", "time_ovarianCA"
)

# Define non-cancer disease event and follow-up time variables.
other_events <- c(
  "alzheimers_disease_group", "copd_group",
  "end_stage_renal_disease_group", "heart_failure_group",
  "hyperthyroidism_group", "hypothyroidism_group",
  "infectious_diseases_group", "parkinsons_disease_group",
  "resp_failure_group", "rheumatoid_arthritis_group",
  "septicaemia_group", "systemic_lupus_erythematosus_group"
)

other_times <- c(
  "survival_time_alzheimers_disease",
  "survival_time_copd",
  "survival_time_end_stage_renal_disease",
  "survival_time_heart_failure",
  "survival_time_hyperthyroidism",
  "survival_time_hypothyroidism",
  "survival_time_infectious_diseases",
  "survival_time_parkinsons_disease",
  "survival_time_resp_failure",
  "survival_time_rheumatoid_arthritis",
  "survival_time_septicaemia",
  "survival_time_systemic_lupus_erythematosus"
)

# Create composite cancer and non-cancer disease variables.
data <- data %>%
  mutate(
    cancer_group = as.integer(
      rowSums(across(all_of(cancer_events)) == 1, na.rm = TRUE) > 0
    ),
    cancer_time = apply(
      across(all_of(cancer_times)),
      1,
      function(x) {
        if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
      }
    )
  ) %>%
  mutate(
    other_group = as.integer(
      rowSums(across(all_of(other_events)) == 1, na.rm = TRUE) > 0
    ),
    other_time = apply(
      across(all_of(other_times)),
      1,
      function(x) {
        if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
      }
    )
  )

# Apply inclusion and exclusion criteria.
data <- data %>%
  filter(BMI_0 > 18.5) %>%
  filter(cancer_time > 0) %>%
  filter(other_time > 0)

write.csv(data, "D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果/raw_data.csv")

# ------------------------------------------------------------------------------
# Analysis 1: Identify significant proteins by follow-up assessment
# ------------------------------------------------------------------------------
categorical_vars <- c("Sex", "Ethnicity", "Edu", "Employed", "Smoke", "Drink")

# Follow-up one: calculate BMI change using BMI_2 and BMI_0.
data <- read.csv("D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果/raw_data.csv")
data <- data %>% mutate(across(all_of(categorical_vars), as.factor))
data$BMIchange <- data$BMI_2 - data$BMI_0

data <- data %>%
  filter(!is.na(BMIchange))

results <- data.frame(
  Variable = character(),
  Coefficient = numeric(),
  P_Value = numeric()
)

for (i in 6:20) {
  formula <- as.formula(
    paste(
      "BMIchange ~",
      colnames(data)[i],
      "+ Age + Sex + Ethnicity + Edu + Employed + Tdi +",
      "Smoke + Drink + Mets + Dietscore + BMI"
    )
  )
  
  fit <- lm(formula, data = data)
  
  coef_value <- summary(fit)$coefficients[2, "Estimate"]
  p_value <- summary(fit)$coefficients[2, "Pr(>|t|)"]
  
  results <- rbind(
    results,
    data.frame(
      Variable = colnames(data)[i],
      Coefficient = coef_value,
      P_Value = p_value
    )
  )
}

# Adjust p-values using the false discovery rate method.
results$FDR_P_Value <- p.adjust(results$P_Value, method = "fdr")

significant_results1 <- subset(results, FDR_P_Value < 0.05)


# Follow-up two: calculate BMI change using BMI_1 and BMI_0.
data <- read.csv("D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果/raw_data.csv")
data <- data %>% mutate(across(all_of(categorical_vars), as.factor))

data$BMIchange <- data$BMI_1 - data$BMI_0

data <- data %>%
  filter(!is.na(BMIchange)) %>%
  filter(is.na(BMI_2))

results <- data.frame(
  Variable = character(),
  Coefficient = numeric(),
  P_Value = numeric()
)

for (i in 6:20) {
  formula <- as.formula(
    paste(
      "BMIchange ~",
      colnames(data)[i],
      "+ Age + Sex + Ethnicity + Edu + Employed + Tdi +",
      "Smoke + Drink + Mets + Dietscore + BMI"
    )
  )
  
  fit <- lm(formula, data = data)
  
  coef_value <- summary(fit)$coefficients[2, "Estimate"]
  p_value <- summary(fit)$coefficients[2, "Pr(>|t|)"]
  
  results <- rbind(
    results,
    data.frame(
      Variable = colnames(data)[i],
      Coefficient = coef_value,
      P_Value = p_value
    )
  )
}

results$FDR_P_Value <- p.adjust(results$P_Value, method = "fdr")

# Select proteins with nominal significance in the validation dataset.
significant_results2 <- subset(results, P_Value < 0.05)

# Identify proteins significant in both follow-up analyses.
common_proteins <- Reduce(
  intersect,
  list(significant_results1$Variable, significant_results2$Variable)
)

# Training set results.
results1_filtered <- subset(
  significant_results1,
  Variable %in% common_proteins
)

# Validation set results.
results2_filtered <- subset(
  significant_results2,
  Variable %in% common_proteins
)

# Merge training and validation results.
final_results <- merge(
  results1_filtered[, c("Variable", "Coefficient", "FDR_P_Value")],
  results2_filtered[, c("Variable", "Coefficient", "P_Value")],
  by = "Variable",
  suffixes = c("_1", "_2")
)

write.csv(final_results, "D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果/weightchange_assessment_time.csv")

# ------------------------------------------------------------------------------
# Analysis 2: Identify significant proteins using a 60:40 train-test split
# ------------------------------------------------------------------------------

data <- read.csv("D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果/raw_data.csv")
data <- data %>% mutate(across(all_of(categorical_vars), as.factor))

# Calculate BMI change using BMI_2 when available; otherwise use BMI_1.
data <- data %>%
  mutate(
    BMIchange = if_else(
      is.na(BMI_2),
      BMI_1 - BMI_0,
      BMI_2 - BMI_0
    )
  )

data <- data %>%
  filter(!is.na(BMIchange))

# Randomly split the data into training and testing datasets.
set.seed(123)

n <- nrow(data)
train_index <- sample(seq_len(n), size = round(0.6 * n))

train_data <- data[train_index, ]
test_data <- data[-train_index, ]


# Training set analysis.
results <- data.frame(
  Variable = character(),
  Coefficient = numeric(),
  P_Value = numeric()
)

for (i in 6:2925) {
  formula <- as.formula(
    paste(
      "BMIchange ~",
      colnames(data)[i],
      "+ Age + Sex + Ethnicity + Edu + Employed + Tdi +",
      "Smoke + Drink + Mets + Dietscore + BMI"
    )
  )
  
  fit <- lm(formula, data = train_data)
  
  coef_value <- summary(fit)$coefficients[2, "Estimate"]
  p_value <- summary(fit)$coefficients[2, "Pr(>|t|)"]
  
  results <- rbind(
    results,
    data.frame(
      Variable = colnames(data)[i],
      Coefficient = coef_value,
      P_Value = p_value
    )
  )
}

results$FDR_P_Value <- p.adjust(results$P_Value, method = "fdr")

significant_results1 <- subset(results, FDR_P_Value < 0.05)


# Testing set analysis.
results <- data.frame(
  Variable = character(),
  Coefficient = numeric(),
  P_Value = numeric()
)

for (i in 6:2925) {
  formula <- as.formula(
    paste(
      "BMIchange ~",
      colnames(data)[i],
      "+ Age + Sex + Ethnicity + Edu + Employed + Tdi +",
      "Smoke + Drink + Mets + Dietscore + BMI"
    )
  )
  
  fit <- lm(formula, data = test_data)
  
  coef_value <- summary(fit)$coefficients[2, "Estimate"]
  p_value <- summary(fit)$coefficients[2, "Pr(>|t|)"]
  
  results <- rbind(
    results,
    data.frame(
      Variable = colnames(data)[i],
      Coefficient = coef_value,
      P_Value = p_value
    )
  )
}

results$FDR_P_Value <- p.adjust(results$P_Value, method = "fdr")

# Select proteins with nominal significance in the testing dataset.
significant_results2 <- subset(results, P_Value < 0.05)

# Identify proteins significant in both training and testing datasets.
common_proteins <- Reduce(
  intersect,
  list(significant_results1$Variable, significant_results2$Variable)
)

# Training set results.
results1_filtered <- subset(
  significant_results1,
  Variable %in% common_proteins
)

# Validation set results.
results2_filtered <- subset(
  significant_results2,
  Variable %in% common_proteins
)

# Merge training and testing results.
final_results <- merge(
  results1_filtered[, c("Variable", "Coefficient", "FDR_P_Value")],
  results2_filtered[, c("Variable", "Coefficient", "P_Value")],
  by = "Variable",
  suffixes = c("_1", "_2")
)

write.csv(final_results, "D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果/weightchange_64分.csv")