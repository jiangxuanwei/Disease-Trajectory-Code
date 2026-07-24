################################################################################
## BMI Loss project — Cox models for incident disease outcomes
##
## Section 1: joint exposure of self-reported weight loss x predicted BMI-loss
##            score tertile (6 mutually exclusive states)
## Section 2: self-reported weight loss alone (binary)
################################################################################

# install.packages("dplyr")
library(dplyr)
library(survival)

## ---------------------------------------------------------------------------
## Paths
## ---------------------------------------------------------------------------

path_weight  <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/UKB_WeightLoss.csv"
path_cov     <- "D:/科研数据/课题数据/UKB分析/UKB数据常用汇总/covariates_impute.csv"
path_outcome <- "D:/科研数据/课题数据/UKB分析/All disease/UBK_disease.csv"
path_score   <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/更新后结果/BMIscore_predict.csv"
dir_out      <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/更新后结果/"

## ---------------------------------------------------------------------------
## Load and merge (done once — both sections use the identical dataset)
## ---------------------------------------------------------------------------

data1     <- read.csv(path_weight)
data3     <- read.csv(path_cov)
outcome   <- read.csv(path_outcome)
datascore <- read.csv(path_score)

# Duplicate IDs would turn each merge into a partial cartesian join and silently
# inflate the sample. Check before merging rather than after.
for (nm in c("data1", "data3", "outcome", "datascore")) {
  d <- get(nm)
  if (any(duplicated(d$n_eid))) {
    stop("Duplicate n_eid found in ", nm, ": ", sum(duplicated(d$n_eid)), " rows.")
  }
}

# Exposure: sign flipped so that a higher value = greater predicted BMI decrease
datascore$BMI_predict_decrease <- datascore$BMI_predict_decrease * (-1)

datacox <- data1 %>%
  merge(data3,     by = "n_eid") %>%
  merge(outcome,   by = "n_eid") %>%
  merge(datascore, by = "n_eid")

cat("Merged analytic file:", nrow(datacox), "participants\n")

# Tertiles of the predicted BMI-loss score (1 = lowest, 3 = highest).
# ntile() assigns NA scores to NA, so tertile boundaries are based on
# non-missing values only.
datacox$score <- ntile(datacox$BMI_predict_decrease, 3)

## Outcome definitions: each event indicator paired with its follow-up time
events1 <- c("alcoholic_liver_disease_group",
             "cirrhosis_group",
             "chronic_kidney_disease_group",
             "diabetes_mellitus_group",
             "event_liverCA",
             "all_cause_death")

times1 <- c("survival_time_alcoholic_liver_disease",
            "survival_time_cirrhosis",
            "survival_time_chronic_kidney_disease",
            "survival_time_diabetes_mellitus",
            "time_liverCA",
            "all_cause_survival_time")

stopifnot(length(events1) == length(times1))
stopifnot(all(c(events1, times1) %in% colnames(datacox)))

covariates <- c("BMI", "Age", "Sex", "Ethnicity", "Edu", "Employed",
                "Tdi", "Smoke", "Drink", "Mets", "Dietscore")

stopifnot(all(covariates %in% colnames(datacox)))

## ---------------------------------------------------------------------------
## Helper: fit one Cox model per outcome and collect the exposure coefficients
## ---------------------------------------------------------------------------

run_cox <- function(data, events, times, exposure_term, covariates) {
  
  res <- lapply(seq_along(events), function(i) {
    event_col <- events[i]
    time_col  <- times[i]
    
    # FIX: the follow-up-time filter is applied to the analysis dataset and the
    # filtered data is what actually reaches coxph(). The original computed
    # `datacox_filtered` from the unrestricted `datacox` and then discarded it —
    # coxph() received the unfiltered `data_loss`, so the filter never ran, and
    # had it run it would have used the wrong sample (no n_2306_0_0 restriction).
    d <- data %>%
      filter(!is.na(.data[[time_col]]),
             .data[[time_col]] > 0,          # t = 0 contributes nothing to the likelihood
             !is.na(.data[[event_col]]))
    
    n_dropped <- nrow(data) - nrow(d)
    if (n_dropped > 0) {
      message(event_col, ": dropped ", n_dropped,
              " rows with missing or non-positive follow-up time.")
    }
    
    if (nrow(d) == 0 || sum(d[[event_col]], na.rm = TRUE) == 0) {
      warning("No events for outcome: ", event_col, " — skipped.")
      return(NULL)
    }
    
    f <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ",
                           exposure_term, " + ",
                           paste(covariates, collapse = " + ")))
    
    fit <- try(coxph(f, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) {
      warning("Cox model failed for outcome: ", event_col)
      return(NULL)
    }
    
    s  <- summary(fit)
    ci <- s$conf.int
    
    # FIX: exposure rows are matched by name instead of the positional indices
    # 1:5 and 1. Positional selection assumes the exposure occupies exactly the
    # first N rows of the coefficient table; if a state cell is empty, or a
    # covariate is character/factor and expands into several rows, the extracted
    # numbers belong to a covariate but get labelled as the exposure.
    keep <- grep(paste0("^", exposure_term), rownames(s$coefficients))
    if (length(keep) == 0) {
      warning("Exposure term not found in model for: ", event_col)
      return(NULL)
    }
    
    data.frame(
      Variable    = rownames(s$coefficients)[keep],
      Coefficient = s$coefficients[keep, "coef"],
      # FIX: take exp(coef) from the summary table rather than recomputing it
      HR          = s$coefficients[keep, "exp(coef)"],
      HR_lower95  = ci[keep, "lower .95"],
      HR_upper95  = ci[keep, "upper .95"],
      p_value     = s$coefficients[keep, "Pr(>|z|)"],
      Event       = event_col,
      N           = fit$n,
      N_events    = fit$nevent,
      stringsAsFactors = FALSE
    )
  })
  
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  
  # FDR across every exposure contrast x outcome combination in this section
  res$FDR_p_value <- p.adjust(res$p_value, method = "fdr")
  res
}


### 1. Self-reported weight loss x predicted BMI-loss score -------------------

# n_2306_0_0: weight change compared with one year ago
#   0 = no change, 1 = gained, 2 = do not know, 3 = lost weight
# Only 0 and 3 are analysed.
cat("\nn_2306_0_0 distribution before restriction:\n")
print(table(datacox$n_2306_0_0, useNA = "ifany"))

data_loss <- datacox %>%
  filter(n_2306_0_0 %in% c(0, 3), !is.na(score))

cat("Section 1 sample:", nrow(data_loss), "participants\n")

# Six mutually exclusive states crossing reported weight change with score
# tertile. Reference = state 1: no reported weight loss AND lowest score tertile.
data_loss <- data_loss %>%
  mutate(state = case_when(
    n_2306_0_0 == 0 & score == 1 ~ 1L,   # reference
    n_2306_0_0 == 0 & score == 2 ~ 2L,
    n_2306_0_0 == 0 & score == 3 ~ 3L,
    n_2306_0_0 == 3 & score == 1 ~ 4L,
    n_2306_0_0 == 3 & score == 2 ~ 5L,
    n_2306_0_0 == 3 & score == 3 ~ 6L,
    TRUE ~ NA_integer_
  ))

# FIX: the factor is built once on the data with an explicit reference level,
# instead of as.factor(state) inline in the formula. This makes the reference
# category visible in the code and gives predictable coefficient row names
# ("state2" ... "state6") for the name-based extraction above.
data_loss$state <- factor(data_loss$state, levels = 1:6)
data_loss <- data_loss %>% filter(!is.na(state))

cat("\nSection 1 state distribution:\n")
print(table(data_loss$state, useNA = "ifany"))

# Empty or very sparse cells make the corresponding HR unstable
if (any(table(data_loss$state) < 10)) {
  warning("At least one state has fewer than 10 participants; ",
          "the corresponding hazard ratios will be unstable.")
}

results_joint <- run_cox(
  data          = data_loss,
  events        = events1,
  times         = times1,
  exposure_term = "state",
  covariates    = covariates
)

write.csv(results_joint, paste0(dir_out, "results_loss.csv"), row.names = FALSE)


### 2. Self-reported weight loss only -----------------------------------------

# The binary exposure does not need the score, so the sample is defined without
# it. This makes Section 2 slightly larger than Section 1.
#
# TO MATCH SECTION 1 EXACTLY, replace the pipeline below with:
#   data_report <- data_loss %>% mutate(reported_loss = ...)
data_report <- datacox %>%
  filter(n_2306_0_0 %in% c(0, 3)) %>%
  # FIX: renamed from `state`. In the original both sections used a variable
  # called `state` with completely different codings (1-6 vs 0/1), which is easy
  # to confuse when the sections are run or edited independently.
  mutate(reported_loss = factor(ifelse(n_2306_0_0 == 3, 1L, 0L),
                                levels = c(0, 1),
                                labels = c("No loss", "Reported loss")))

cat("\nSection 2 sample:", nrow(data_report), "participants\n")
print(table(data_report$reported_loss, useNA = "ifany"))

results_report <- run_cox(
  data          = data_report,
  events        = events1,
  times         = times1,
  exposure_term = "reported_loss",
  covariates    = covariates
)

write.csv(results_report, paste0(dir_out, "results_report_loss.csv"), row.names = FALSE)


### Optional: proportional hazards check --------------------------------------
## Cox models assume proportional hazards. Worth checking at least for the
## primary outcome before reporting. A significant global test means the HR is
## not constant over follow-up and a time-varying term or stratification may be
## needed.
# d_check <- data_loss %>% filter(all_cause_survival_time > 0)
# fit_check <- coxph(
#   Surv(all_cause_survival_time, all_cause_death) ~ state + BMI + Age + Sex +
#     Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore,
#   data = d_check
# )
# print(cox.zph(fit_check))


### Optional: global FDR across both sections ---------------------------------
## Each section corrects within itself (30 tests in Section 1, 6 in Section 2).
## For a single study-wide correction:
# all_results <- bind_rows(
#   results_joint  %>% mutate(Model = "Joint"),
#   results_report %>% mutate(Model = "Reported loss only")
# )
# all_results$FDR_p_global <- p.adjust(all_results$p_value, method = "fdr")
# write.csv(all_results, paste0(dir_out, "results_all.csv"), row.names = FALSE)