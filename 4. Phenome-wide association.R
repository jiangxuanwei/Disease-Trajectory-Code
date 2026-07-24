### Proteomics-predicted BMI decrease and incident disease risk (UKB)
### Cox proportional hazards models across multiple outcomes

library(dplyr)
library(survival)

# ---- Paths ----------------------------------------------------------------
base_dir  <- "D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果"
raw_file  <- file.path(base_dir, "raw_data.csv")
score_file <- file.path(base_dir, "BMIscore_predict.csv")
out_file  <- file.path(base_dir, "proteomics_BMI decrease per sd.csv")

# ---- Load and merge data --------------------------------------------------
data      <- read.csv(raw_file)
datascore <- read.csv(score_file)
datacox   <- merge(data, datascore, by = "n_eid")

# ---- Exposure -------------------------------------------------------------
# Reverse sign so that higher values indicate greater predicted BMI decrease
datacox$BMI_predict_decrease <- datacox$BMI_predict_decrease * (-1)

# Standardized exposure (per 1 SD increase)
datacox$score <- datacox$BMI_predict_decrease / sd(datacox$BMI_predict_decrease)

# ---- Follow-up time -------------------------------------------------------
# For participants without a disease-specific survival time, use death follow-up time
datacox <- datacox %>%
  mutate(across(starts_with("survival_time_"), ~ ifelse(is.na(.), death2_time, .)))

# ---- Outcome definitions --------------------------------------------------
events <- c("cvd_group", "cad_group", "hypertension_group", "heart_failure_group", "arrhythmia_group",
            "stroke_group", "aaa_group", "cardiomyopathy_group", "pericardial_group", "nrvd_group",
            "pe_group", "pulm_hd_group", "pvd_group", "mvd_group", "phlebitis_group", "All_CVD_events_group",
            "asthma_group", "copd_group", "pneumonia_group", "bronchiectasis_group", "resp_failure_group",
            "other_resp_disorders_group", "all_resp_diseases_group", "liver_disease_group", "alcoholic_liver_disease_group",
            "chronic_hepatitis_group", "cirrhosis_group", "fatty_liver_group", "gastric_ulcer_group",
            "peptic_ulcer_group", "gastro_oesophageal_reflux_disease_group", "pancreatitis_group",
            "gallbladder_biliary_tract_diseases_group", "ibs_group", "digestive_system_diseases_group",
            "acute_kidney_injury_group", "chronic_kidney_disease_group", "end_stage_renal_disease_group",
            "kidney_stones_group", "prostatitis_group", "urinary_system_diseases_group",
            "alzheimers_disease_group", "parkinsons_disease_group", "epilepsy_group",
            "multiple_sclerosis_group", "migraine_group", "autonomic_nervous_system_disorders_group",
            "facial_nerve_disorders_group", "neuro_system_diseases_group", "hyperlipidemia_group",
            "diabetes_mellitus_group", "hyperthyroidism_group", "hypothyroidism_group",
            "iodine_deficiency_thyroid_disease_group", "obesity_group", "endocrine_metabolic_diseases_group",
            "osteoporosis_group", "osteoarthritis_group", "musculoskeletal_diseases_group",
            "immunodeficiency_group", "rheumatoid_arthritis_group", "systemic_lupus_erythematosus_group",
            "immunological_diseases_group", "bipolar_disorder_group", "depression_group",
            "anxiety_disorders_group", "schizophrenia_group", "adhd_group", "alcohol_problems_group",
            "mental_health_disorders_group", "hearing_loss_group", "visual_impairment_blindness_group",
            "cataract_group", "glaucoma_group", "sensory_system_disorders_group", "anemia_group",
            "thrombocytopaenia_group", "thrombophilia_group", "blood_disorders_group", "septicaemia_group",
            "skin_infection_group", "rheumatic_fever_group", "infectious_diseases_group",
            "actinic_keratosis_group", "dermatitis_group", "gout_group", "psoriasis_group",
            "urticaria_group", "lichen_planus_group", "skin_diseases_group", "all_cause_death",
            "event_LungCA", "event_BreastCA", "event_ColorectalCA",
            "event_ProstaticCA", "event_GastroCA", "event_liverCA", "event_EsophagealCA",
            "event_PancreaticCA", "event_CervicalCA", "event_ovarianCA")

times <- c("survival_time_cvd", "survival_time_cad", "survival_time_hypertension", "survival_time_heart_failure",
           "survival_time_arrhythmia", "survival_time_stroke", "survival_time_aaa", "survival_time_cardiomyopathy",
           "survival_time_pericardial", "survival_time_nrvd", "survival_time_pe", "survival_time_pulm_hd",
           "survival_time_pvd", "survival_time_mvd", "survival_time_phlebitis", "survival_time_All_CVD_events",
           "survival_time_asthma", "survival_time_copd", "survival_time_pneumonia", "survival_time_bronchiectasis",
           "survival_time_resp_failure", "survival_time_other_resp_disorders", "survival_time_all_resp_diseases",
           "survival_time_liver_disease", "survival_time_alcoholic_liver_disease", "survival_time_chronic_hepatitis",
           "survival_time_cirrhosis", "survival_time_fatty_liver", "survival_time_gastric_ulcer",
           "survival_time_peptic_ulcer", "survival_time_gastro_oesophageal_reflux_disease", "survival_time_pancreatitis",
           "survival_time_gallbladder_biliary_tract_diseases", "survival_time_ibs", "survival_time_digestive_system_diseases",
           "survival_time_acute_kidney_injury", "survival_time_chronic_kidney_disease",
           "survival_time_end_stage_renal_disease", "survival_time_kidney_stones", "survival_time_prostatitis",
           "survival_time_urinary_system_diseases", "survival_time_alzheimers_disease",
           "survival_time_parkinsons_disease", "survival_time_epilepsy", "survival_time_multiple_sclerosis",
           "survival_time_migraine", "survival_time_autonomic_nervous_system_disorders",
           "survival_time_facial_nerve_disorders", "survival_time_neuro_system_diseases", "survival_time_hyperlipidemia",
           "survival_time_diabetes_mellitus", "survival_time_hyperthyroidism", "survival_time_hypothyroidism",
           "survival_time_iodine_deficiency_thyroid_disease", "survival_time_obesity",
           "survival_time_endocrine_metabolic_diseases", "survival_time_osteoporosis",
           "survival_time_osteoarthritis", "survival_time_musculoskeletal_diseases",
           "survival_time_immunodeficiency", "survival_time_rheumatoid_arthritis",
           "survival_time_systemic_lupus_erythematosus", "survival_time_immunological_diseases",
           "survival_time_bipolar_disorder", "survival_time_depression", "survival_time_anxiety_disorders",
           "survival_time_schizophrenia", "survival_time_adhd", "survival_time_alcohol_problems",
           "survival_time_mental_health_disorders", "survival_time_hearing_loss",
           "survival_time_visual_impairment_blindness", "survival_time_cataract", "survival_time_glaucoma",
           "survival_time_sensory_system_disorders", "survival_time_anemia", "survival_time_thrombocytopaenia",
           "survival_time_thrombophilia", "survival_time_blood_disorders", "survival_time_septicaemia",
           "survival_time_skin_infection", "survival_time_rheumatic_fever", "survival_time_infectious_diseases",
           "survival_time_actinic_keratosis", "survival_time_dermatitis", "survival_time_gout",
           "survival_time_psoriasis", "survival_time_urticaria", "survival_time_lichen_planus",
           "survival_time_skin_diseases", "all_cause_survival_time",
           "time_LungCA", "time_BreastCA", "time_ColorectalCA", "time_ProstaticCA", "time_GastroCA",
           "time_liverCA", "time_EsophagealCA", "time_PancreaticCA", "time_CervicalCA", "time_ovarianCA")

stopifnot(length(events) == length(times))

# ---- Covariates -----------------------------------------------------------
covariates <- c("BMI_0", "Age", "Sex", "Ethnicity", "Edu", "Employed",
                "Tdi", "Smoke", "Drink", "Mets", "Dietscore")

# ---- Cox models across all outcomes ---------------------------------------
results <- data.frame()

for (i in seq_along(events)) {
  event_col <- events[i]
  time_col  <- times[i]
  
  # Keep participants with non-negative follow-up time
  datacox_filtered <- datacox %>% filter(.data[[time_col]] >= 0)
  
  formula <- as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ BMI_predict_decrease + ",
           paste(covariates, collapse = " + "))
  )
  
  cox_model   <- coxph(formula, data = datacox_filtered)
  summary_cox <- summary(cox_model)
  
  result_df <- data.frame(
    Event        = event_col,
    Variable     = rownames(summary_cox$coefficients)[1],
    N_total      = summary_cox$n,            # participants included in the model (complete cases)
    N_events     = summary_cox$nevent,       # incident cases
    N_nonevents  = summary_cox$n - summary_cox$nevent,
    N_filtered   = nrow(datacox_filtered),   # rows after time filter, before listwise deletion
    N_missing    = nrow(datacox_filtered) - summary_cox$n,
    Coefficient  = summary_cox$coefficients[1, "coef"],
    HR           = exp(summary_cox$coefficients[1, "coef"]),
    HR_lower95   = summary_cox$conf.int[1, "lower .95"],
    HR_upper95   = summary_cox$conf.int[1, "upper .95"],
    p_value      = summary_cox$coefficients[1, "Pr(>|z|)"]
  )
  
  results <- rbind(results, result_df)
}

# ---- Multiple testing correction ------------------------------------------
results$FDR_p_value <- p.adjust(results$p_value, method = "fdr")

# ---- Export ---------------------------------------------------------------
write.csv(results, out_file, row.names = FALSE)