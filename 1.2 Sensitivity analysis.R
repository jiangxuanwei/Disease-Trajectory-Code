library(dplyr)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

raw_data_file <- "D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果/raw_data.csv"
med_data_file <- "D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果/weight_loss_med.csv"
bmi_file <- "D:/科研数据/课题数据/1. UKB数据课题/UKB数据常用汇总/Obesity data/BMI.csv"
protein_file <- "D:/科研数据/课题数据/1. UKB数据课题/UKB数据常用汇总/Omics data/protein_UKB_filled.csv"
covariate_file <- "D:/科研数据/课题数据/1. UKB数据课题/UKB数据常用汇总/covariates_impute.csv"
disease_file <- "D:/科研数据/课题数据/1. UKB数据课题/All disease/UBK_disease.csv"

out_dir <- "D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/投稿/sub for Metabolism/第一轮重投/补充结果"

# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

proteins <- c("ache", "adamts8", "angptl7", "cd93", "comp", "dlk1", "fap", "gcg",
              "gfra3", "igsf3", "ism1", "itgb6", "prtg", "ret", "thbs4")

covariates_main <- c("Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi",
                     "Smoke", "Drink", "Mets", "Dietscore", "BMI")

categorical_vars <- c("Sex", "Ethnicity", "Edu", "Employed", "Smoke", "Drink")

# ------------------------------------------------------------------------------
# Function: protein-BMIchange association
# ------------------------------------------------------------------------------

run_protein_bmichange_association <- function(input_data, output_suffix) {
  input_data <- input_data %>% mutate(across(any_of(categorical_vars), as.factor))
  results <- data.frame(Variable = proteins, N = NA, Coefficient = NA, SE = NA, P_Value = NA)
  
  for (i in seq_along(proteins)) {
    protein <- proteins[i]
    required_vars <- c("BMIchange", protein, covariates_main)
    
    if (all(required_vars %in% colnames(input_data))) {
      model_data <- input_data %>% select(all_of(required_vars)) %>% na.omit()
      results$N[i] <- nrow(model_data)
      
      if (nrow(model_data) > 30 && length(unique(model_data[[protein]])) > 1) {
        formula <- as.formula(
          paste("BMIchange ~", protein,
                "+ Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore + BMI")
        )
        
        fit <- tryCatch(lm(formula, data = model_data), error = function(e) NULL)
        
        if (!is.null(fit)) {
          coef_table <- summary(fit)$coefficients
          if (protein %in% rownames(coef_table)) {
            results$Coefficient[i] <- coef_table[protein, "Estimate"]
            results$SE[i] <- coef_table[protein, "Std. Error"]
            results$P_Value[i] <- coef_table[protein, "Pr(>|t|)"]
          }
        }
      }
    }
  }
  
  results$FDR_P_Value <- p.adjust(results$P_Value, method = "fdr")
  write.csv(results, file.path(out_dir, paste0("protein_bmichange_", output_suffix, ".csv")), row.names = FALSE)
  return(results)
}

# ------------------------------------------------------------------------------
# 1. Sensitivity analysis using raw_data: lifestyle and medication
# ------------------------------------------------------------------------------

raw_data <- read.csv(raw_data_file) %>%
  mutate(BMIchange = if_else(is.na(BMI_2), BMI_1 - BMI_0, BMI_2 - BMI_0)) %>%
  filter(!is.na(BMIchange))

# ---- Merge weight loss medication data ----
med_data <- read.csv(med_data_file)

id_var <- intersect(names(raw_data), names(med_data))[1]

med_data <- med_data %>%
  select(all_of(c(id_var, "weight_loss_med"))) %>%
  distinct(.data[[id_var]], .keep_all = TRUE)

raw_data <- raw_data %>%
  left_join(med_data, by = id_var)

raw_data <- raw_data %>%
  mutate(
    smoke_ideal = case_when(Smoke == 0 ~ 1, Smoke > 0 ~ 0, TRUE ~ NA_real_),
    pa_ideal = case_when(Mets >= 600 ~ 1, Mets < 600 ~ 0, TRUE ~ NA_real_),
    diet_ideal = case_when(Dietscore >= 4 ~ 1, Dietscore < 4 ~ 0, TRUE ~ NA_real_),
    sleep_ideal = case_when(Sleeptime >= 7 & Sleeptime <= 9 ~ 1,
                            Sleeptime < 7 | Sleeptime > 9 ~ 0,
                            TRUE ~ NA_real_),
    lifestyle_score = case_when(
      !is.na(smoke_ideal) & !is.na(pa_ideal) & !is.na(diet_ideal) & !is.na(sleep_ideal) ~
        smoke_ideal + pa_ideal + diet_ideal + sleep_ideal,
      TRUE ~ NA_real_
    )
  )

results_original <- run_protein_bmichange_association(raw_data, "original_raw_data")

data_exclude_lifestyle_score4 <- raw_data %>%
  filter(is.na(lifestyle_score) | lifestyle_score != 4)

results_exclude_lifestyle <- run_protein_bmichange_association(
  data_exclude_lifestyle_score4,
  "exclude_lifestyle_score4"
)

# ---- Exclude metformin users AND weight loss medication users ----
data_exclude_med <- raw_data %>%
  filter((is.na(med_Metformin) | med_Metformin != 1) &
           (is.na(weight_loss_med) | weight_loss_med != 1))

results_exclude_med <- run_protein_bmichange_association(
  data_exclude_med,
  "exclude_metformin_and_weightlossmed"
)

# ------------------------------------------------------------------------------
# 2. Sensitivity analysis using original disease files: exclude diseases within 2 years
# ------------------------------------------------------------------------------

data1 <- read.csv(bmi_file)
data2 <- read.csv(protein_file)
data3 <- read.csv(covariate_file)
track <- read.csv(disease_file)

# If raw_data used standardized proteins, keep this line active.
# data2[, 2:2921] <- scale(log2(data2[, 2:2921] + 20))

disease_data <- merge(data1, data2, by = "n_eid")
disease_data <- merge(disease_data, data3, by = "n_eid")
disease_data <- merge(disease_data, track, by = "n_eid")

disease_data <- disease_data %>%
  mutate(across(starts_with("survival_time_"), ~ ifelse(is.na(.), death2_time, .)))

cancer_events <- c("event_LungCA", "event_BreastCA", "event_ColorectalCA",
                   "event_ProstaticCA", "event_GastroCA", "event_liverCA",
                   "event_EsophagealCA", "event_PancreaticCA",
                   "event_CervicalCA", "event_ovarianCA")

cancer_times <- c("time_LungCA", "time_BreastCA", "time_ColorectalCA",
                  "time_ProstaticCA", "time_GastroCA", "time_liverCA",
                  "time_EsophagealCA", "time_PancreaticCA",
                  "time_CervicalCA", "time_ovarianCA")

other_events <- c("alzheimers_disease_group", "copd_group", "end_stage_renal_disease_group",
                  "heart_failure_group", "hyperthyroidism_group", "hypothyroidism_group",
                  "infectious_diseases_group", "parkinsons_disease_group", "resp_failure_group",
                  "rheumatoid_arthritis_group", "septicaemia_group", "systemic_lupus_erythematosus_group")

other_times <- c("survival_time_alzheimers_disease", "survival_time_copd",
                 "survival_time_end_stage_renal_disease", "survival_time_heart_failure",
                 "survival_time_hyperthyroidism", "survival_time_hypothyroidism",
                 "survival_time_infectious_diseases", "survival_time_parkinsons_disease",
                 "survival_time_resp_failure", "survival_time_rheumatoid_arthritis",
                 "survival_time_septicaemia", "survival_time_systemic_lupus_erythematosus")

disease_data <- disease_data %>%
  mutate(
    cancer_group = as.integer(rowSums(across(all_of(cancer_events)) == 1, na.rm = TRUE) > 0),
    cancer_time = apply(across(all_of(cancer_times)), 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)),
    other_group = as.integer(rowSums(across(all_of(other_events)) == 1, na.rm = TRUE) > 0),
    other_time = apply(across(all_of(other_times)), 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
  ) %>%
  filter(BMI_0 > 18.5) %>%
  filter(cancer_time > 0) %>%
  filter(other_time > 0) %>%
  mutate(BMIchange = if_else(is.na(BMI_2), BMI_1 - BMI_0, BMI_2 - BMI_0)) %>%
  filter(!is.na(BMIchange))

data_exclude_cancer_2y <- disease_data %>%
  filter(is.na(cancer_time) | cancer_time > 2)

results_exclude_cancer_2y <- run_protein_bmichange_association(
  data_exclude_cancer_2y,
  "exclude_cancer_within_2y"
)

data_exclude_other_disease_2y <- disease_data %>%
  filter(is.na(other_time) | other_time > 2)

results_exclude_other_disease_2y <- run_protein_bmichange_association(
  data_exclude_other_disease_2y,
  "exclude_other_disease_within_2y"
)

data_exclude_any_disease_2y <- disease_data %>%
  filter((is.na(cancer_time) | cancer_time > 2) &
           (is.na(other_time) | other_time > 2))

results_exclude_any_disease_2y <- run_protein_bmichange_association(
  data_exclude_any_disease_2y,
  "exclude_any_disease_within_2y"
)

# ------------------------------------------------------------------------------
# 3. Sensitivity analysis using original disease files: exclude diseases within 5 years
# ------------------------------------------------------------------------------

data_exclude_cancer_5y <- disease_data %>%
  filter(is.na(cancer_time) | cancer_time > 5)

results_exclude_cancer_5y <- run_protein_bmichange_association(
  data_exclude_cancer_5y,
  "exclude_cancer_within_5y"
)

data_exclude_other_disease_5y <- disease_data %>%
  filter(is.na(other_time) | other_time > 5)

results_exclude_other_disease_5y <- run_protein_bmichange_association(
  data_exclude_other_disease_5y,
  "exclude_other_disease_within_5y"
)

data_exclude_any_disease_5y <- disease_data %>%
  filter((is.na(cancer_time) | cancer_time > 5) &
           (is.na(other_time) | other_time > 5))

results_exclude_any_disease_5y <- run_protein_bmichange_association(
  data_exclude_any_disease_5y,
  "exclude_any_disease_within_5y"
)