library(dplyr)

# =========================================================
# 1. Read and prepare data
# =========================================================

data <- read.csv("D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/更新后结果/raw_data.csv")
data1 <- read.csv("D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/UKB_WeightLoss.csv")

data <- merge(data, data1, by = "n_eid") %>%
  mutate(BMIchange = if_else(is.na(BMI_2), BMI_1 - BMI_0, BMI_2 - BMI_0)) %>%
  filter(!is.na(BMIchange))

# Convert categorical covariates to factors.
categorical_vars <- c("Sex", "Ethnicity", "Edu", "Employed", "Smoke", "Drink")
data <- data %>% mutate(across(all_of(categorical_vars), as.factor))

proteins <- c("ache", "adamts8", "angptl7", "cd93", "comp", "dlk1", "fap", "gcg",
              "gfra3", "igsf3", "ism1", "itgb6", "prtg", "ret", "thbs4")

covariates_main <- c("Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi",
                     "Smoke", "Drink", "Mets", "Dietscore", "BMI")

out_dir <- "D:/科研数据/课题数据/1. UKB数据课题/BMI Loss 课题/投稿/sub for Metabolism/第一轮重投/补充结果"


# =========================================================
# 2. Function: subgroup analysis
# =========================================================

run_subgroup_analysis <- function(input_data, output_suffix) {
  
  age_cutoffs <- c(50, 55, 60, 65)
  
  age_subgroups <- unlist(
    lapply(age_cutoffs, function(cutoff) {
      setNames(
        list(input_data %>% filter(Age < cutoff),
             input_data %>% filter(Age >= cutoff)),
        c(paste0("Age_lt", cutoff), paste0("Age_ge", cutoff))
      )
    }),
    recursive = FALSE
  )
  
  other_subgroups <- list(
    BMI_lt25 = input_data %>% filter(BMI < 25),
    BMI_25_30 = input_data %>% filter(BMI >= 25 & BMI < 30),
    BMI_ge30 = input_data %>% filter(BMI >= 30),
    weight_stable = input_data %>% filter(n_2306_0_0 == 0),
    weight_loss = input_data %>% filter(n_2306_0_0 == 3)
  )
  
  subgroups <- c(age_subgroups, other_subgroups)
  all_results <- list()
  
  for (group_name in names(subgroups)) {
    sub_data <- subgroups[[group_name]]
    
    results <- data.frame(Variable = proteins, N = NA, Coefficient = NA, SE = NA, P_Value = NA)
    
    for (i in seq_along(proteins)) {
      protein <- proteins[i]
      required_vars <- c("BMIchange", protein, covariates_main)
      
      if (all(required_vars %in% colnames(sub_data))) {
        model_data <- sub_data %>% select(all_of(required_vars)) %>% na.omit()
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
    colnames(results)[2:6] <- paste0(colnames(results)[2:6], "_", group_name)
    all_results[[group_name]] <- results
  }
  
  final_results <- Reduce(function(x, y) merge(x, y, by = "Variable", all = TRUE), all_results)
  
  write.csv(
    final_results,
    file.path(out_dir, paste0("subgroup_analysis_", output_suffix, ".csv")),
    row.names = FALSE
  )
  
  return(final_results)
}


# =========================================================
# 3. Run subgroup analysis
# =========================================================

results_subgroup <- run_subgroup_analysis(input_data = data, output_suffix = "original")
