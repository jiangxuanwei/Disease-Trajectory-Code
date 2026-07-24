################################################################################
## BMI Loss project — discrimination of nested prediction models
##
## Compares four nested Cox models per outcome using time-dependent AUC at
## 15 years:
##   fs1 Base    : age, sex, BMI
##   fs2 BMI     : base + lifestyle factors
##   fs3 Panel   : lifestyle + clinical biochemistry panel
##   fs4 Prot    : panel + 15 selected proteins
################################################################################

rm(list = ls())

library(riskRegression)
library(survival)
library(prodlim)
library(pec)
library(dplyr)

## ---------------------------------------------------------------------------
## Paths
## ---------------------------------------------------------------------------

dir_base <- "D:/科研数据/课题数据/1. UKB数据课题/"
dir_out  <- paste0(dir_base, "BMI Loss 课题/更新后结果/")

path_cov     <- paste0(dir_base, "UKB数据常用汇总/covariates_impute.csv")
path_outcome <- paste0(dir_base, "All disease/UBK_disease.csv")
path_protein <- paste0(dir_base, "UKB数据常用汇总/Omics data/protein_UKB_filled.csv")
path_panel   <- paste0(dir_base, "BMI Loss 课题/PANEL.csv")
path_bmi     <- paste0(dir_base, "UKB数据常用汇总/Obesity data/BMI.csv")
path_score   <- paste0(dir_base, "BMI Loss 课题/更新后结果/BMIscore_predict.csv")

## ---------------------------------------------------------------------------
## Load data
## ---------------------------------------------------------------------------

data3     <- read.csv(path_cov)          # covariates (imputed)
outcome   <- read.csv(path_outcome)      # disease events and follow-up times
data2     <- read.csv(path_protein)      # Olink protein panel
panel     <- read.csv(path_panel)        # clinical biochemistry panel
BMI       <- read.csv(path_bmi)          # repeated BMI measurements
datascore <- read.csv(path_score)        # predicted BMI-loss score

# Protein columns: log2-transform with an offset of 20 to keep values positive,
# then z-standardise so all proteins share a common scale.
# NOTE: the hard-coded range 2:2921 breaks silently if the file's column count
# ever changes. Selecting by position of the ID column is safer.
protein_cols <- setdiff(colnames(data2), "n_eid")
#data2[protein_cols] <- scale(log2(data2[protein_cols] + 20))

## ---------------------------------------------------------------------------
## Merge into the analytic dataset
## ---------------------------------------------------------------------------

datacox <- data3 %>%
  merge(outcome,   by = "n_eid") %>%
  merge(data2,     by = "n_eid") %>%
  merge(panel,     by = "n_eid") %>%
  merge(BMI,       by = "n_eid") %>%
  merge(datascore, by = "n_eid")

cat("After merging:", nrow(datacox), "participants\n")

# Restrict to participants without repeat BMI measurements at instances 1 and 2.
# See the note below the script: as written this almost certainly empties the
# data. Verify what BMI_1 / BMI_2 encode before running.
datacox <- datacox %>% filter(is.na(BMI_1) & is.na(BMI_2))

cat("After BMI restriction:", nrow(datacox), "participants\n")
stopifnot(nrow(datacox) > 0)

## ---------------------------------------------------------------------------
## Outcome and model definitions
## ---------------------------------------------------------------------------

events1 <- c("alcoholic_liver_disease_group",
             "cirrhosis_group",
             "chronic_kidney_disease_group",
             "diabetes_mellitus_group",
             "event_liverCA")

times1 <- c("survival_time_alcoholic_liver_disease",
            "survival_time_cirrhosis",
            "survival_time_chronic_kidney_disease",
            "survival_time_diabetes_mellitus",
            "time_liverCA")

# The .x suffixes come from duplicated column names across the merged files.
# Renaming those columns after merging would make these formulas far easier to
# read and remove the risk of silently picking up the wrong copy.
Lifestyle <- "Age.x + Sex.x + BMI.x + Smoke.x + Drink.x + Sleeptime + Dietscore.x"

Lifestyle_add <- paste(Lifestyle,
                       "+ TC + HDL + glucose + HbA1c + AST + ALT + AP + Albumin",
                       "+ Creatinine + Cystatine_C + Urea + Urate")

# Fifteen proteins selected for the final model
proteins <- c("ache", "adamts8", "angptl7", "cd93", "comp", "dlk1", "fap",
              "gcg", "gfra3", "igsf3", "ism1", "itgb6", "prtg", "ret", "thbs4")

Protein_model <- paste(Lifestyle_add, "+", paste(proteins, collapse = " + "))

# Prediction horizon for the time-dependent AUC
horizon <- 15

## ---------------------------------------------------------------------------
## Fit and compare models for each outcome
## ---------------------------------------------------------------------------

resultAUC <- list()
resultp   <- list()

for (i in seq_along(events1)) {
  
  outcome_var <- events1[i]
  time_var    <- times1[i]
  
  cat("\n=== Outcome:", outcome_var, "===\n")
  
  # Complete-case subset shared by all four models. Score() compares models on
  # a common dataset, so fitting each model on its own complete cases would make
  # the AUC contrasts incomparable.
  model_vars <- unique(c(time_var, outcome_var, "Age.x", "Sex.x", "BMI.x",
                         "Smoke.x", "Drink.x", "Sleeptime", "Dietscore.x",
                         "TC", "HDL", "glucose", "HbA1c", "AST", "ALT", "AP",
                         "Albumin", "Creatinine", "Cystatine_C", "Urea",
                         "Urate", proteins))
  
  d <- datacox %>%
    filter(if_all(all_of(model_vars), ~ !is.na(.))) %>%
    filter(.data[[time_var]] > 0)
  
  cat("Analytic sample:", nrow(d),
      "| events:", sum(d[[outcome_var]], na.rm = TRUE), "\n")
  
  # Time-dependent AUC at 15 years needs adequate follow-up and events at that
  # horizon; otherwise the estimate is unstable or undefined.
  if (max(d[[time_var]], na.rm = TRUE) < horizon) {
    warning("Maximum follow-up is shorter than the ", horizon,
            "-year horizon for ", outcome_var, "; skipped.")
    next
  }
  if (sum(d[[outcome_var]] == 1 & d[[time_var]] <= horizon, na.rm = TRUE) < 20) {
    warning("Fewer than 20 events before ", horizon, " years for ",
            outcome_var, "; AUC estimates will be unstable.")
  }
  
  surv_lhs <- paste0("Surv(", time_var, ", ", outcome_var, ") ~ ")
  
  # x = TRUE stores the design matrix, which Score() requires
  fs1 <- coxph(as.formula(paste0(surv_lhs, "Age.x + Sex.x + BMI.x")),
               data = d, x = TRUE)
  fs2 <- coxph(as.formula(paste0(surv_lhs, Lifestyle)),
               data = d, x = TRUE)
  fs3 <- coxph(as.formula(paste0(surv_lhs, Lifestyle_add)),
               data = d, x = TRUE)
  fs4 <- coxph(as.formula(paste0(surv_lhs, Protein_model)),
               data = d, x = TRUE)
  
  # Score() evaluates discrimination against a null (marginal) model and
  # produces pairwise contrasts between the listed models
  score_formula <- as.formula(paste0("Hist(", time_var, ", ", outcome_var, ") ~ 1"))
  
  xs <- try(
    Score(list(Base_model = fs1,
               BMI_model  = fs2,
               Panel      = fs3,
               Prot       = fs4),
          formula = score_formula,
          data    = d,
          times   = horizon,
          plots   = "roc",
          metrics = "auc"),
    silent = TRUE
  )
  
  if (inherits(xs, "try-error")) {
    warning("Score() failed for outcome: ", outcome_var)
    next
  }
  
  # Tag results with the outcome so the combined output stays interpretable
  auc_tab <- as.data.frame(xs$AUC$score)
  auc_tab$Outcome <- outcome_var
  resultAUC[[outcome_var]] <- auc_tab
  
  contrast_tab <- as.data.frame(xs$AUC$contrasts)
  contrast_tab$Outcome <- outcome_var
  resultp[[outcome_var]] <- contrast_tab
  
  # ROC curves; uncomment the pdf()/dev.off() pair to write them to file
  # pdf(paste0(dir_out, "ROC_", outcome_var, ".pdf"), onefile = FALSE)
  plotROC(xs, col = c("#7F7F80", "#55719D", "#EF7C1C", "#BF1E20"))
  # dev.off()
}

## ---------------------------------------------------------------------------
## Export
## ---------------------------------------------------------------------------

# FIX: the original passed a list of data frames straight to write.csv(), which
# coerces them into a single wide row-recycled table with mangled column names.
# bind_rows() produces a proper long table with one row per model per outcome.
auc_all      <- bind_rows(resultAUC)
contrast_all <- bind_rows(resultp)

# Multiple-testing correction across all model contrasts and outcomes
if ("p" %in% colnames(contrast_all)) {
  contrast_all$FDR_p <- p.adjust(contrast_all$p, method = "fdr")
}

write.csv(auc_all,      paste0(dir_out, "AUC 15 y.csv"),   row.names = FALSE)
write.csv(contrast_all, paste0(dir_out, "AUC-p 15 y.csv"), row.names = FALSE)