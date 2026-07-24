################################################################################
## BMI Loss project — associations between predicted BMI decrease and
## multi-organ phenotypes (brain, cerebellum, vascular, biochemistry,
## body composition, cardiac)
##
## Revised version: fixes to scaled-variable selection, log transforms,
## left/right combination, coefficient extraction, and factor handling.
################################################################################

library(dplyr)

## ---------------------------------------------------------------------------
## Shared paths and helpers
## ---------------------------------------------------------------------------

path_score  <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/更新后结果/BMIscore_predict.csv"
path_cov    <- "D:/科研数据/课题数据/UKB分析/UKB数据常用汇总/covariates_impute.csv"
path_center <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/center.csv"
dir_raw     <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/原始数据/"
dir_out     <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/更新后结果/"

# Exposure: higher value = greater BMI decrease
load_score <- function() {
  d <- read.csv(path_score)
  d$BMI_predict_decrease <- d$BMI_predict_decrease * (-1)
  d
}

# log2-transform then z-standardise, returning a plain vector (not a 1-col matrix).
# Non-positive values become NA rather than -Inf.
log2_z <- function(x) {
  x <- as.numeric(x)
  x[!is.na(x) & x <= 0] <- NA
  as.vector(scale(log2(x)))
}

# Row sums that return NA when any input is missing (avoids treating missing as 0)
rowsum_complete <- function(df) {
  out <- rowSums(df)          # na.rm = FALSE: any NA propagates
  as.numeric(out)
}

# Extract the exposure coefficient by NAME, not by row position
extract_exposure <- function(fit, term = "BMI_predict_decrease", label) {
  cs <- summary(fit)$coefficients
  if (!term %in% rownames(cs)) return(NULL)
  data.frame(
    Variable = label,
    Beta     = cs[term, "Estimate"],
    SE       = cs[term, "Std. Error"],
    P_value  = cs[term, "Pr(>|t|)"],
    N        = length(fit$residuals),
    stringsAsFactors = FALSE
  )
}

# Fit one model per outcome and collect exposure coefficients
run_models <- function(data, y_vars, predictors) {
  res <- lapply(y_vars, function(y_var) {
    f <- as.formula(paste(y_var, "~", paste(predictors, collapse = " + ")))
    fit <- try(lm(f, data = data, na.action = na.omit), silent = TRUE)
    if (inherits(fit, "try-error")) {
      warning("Model failed for outcome: ", y_var)
      return(NULL)
    }
    extract_exposure(fit, label = y_var)
  })
  res <- do.call(rbind, res)
  res$FDR_P <- p.adjust(res$P_value, method = "fdr")
  res
}

covariates_base <- c("BMI_predict_decrease", "BMI", "Age", "Sex", "Ethnicity",
                     "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore")

# Covariate table used by several sections — loaded once, up front
data3 <- read.csv(path_cov)


### 1. Brain MRI (cerebral regions) -------------------------------------------

datascore <- load_score()
brain     <- read.csv(paste0(dir_raw, "MRI-out.csv"))
center    <- read.csv(path_center)

# FIX: dropped the blanket na.omit(brain). Row-wise deletion across the whole
# IDP table removed participants missing ANY imaging variable, which shrank the
# sample non-randomly. lm(na.action = na.omit) already handles missingness
# per-model, using the maximum available sample for each outcome.

data <- merge(datascore, brain,  by = "n_eid")
data <- merge(data,      center, by = "n_eid")
data <- merge(data,      data3,  by = "n_eid")

head_scaling <- "X25000.2.0"   # T1 head-size scaling factor

## Group 1: raw volumes, multiplied by the head-size scaling factor
group1_ids   <- c("26526", "25003", "25892")
group1_cols  <- paste0("X", group1_ids, ".2.0")
group1_names <- c("brain_stem",
                  "ventricular_cerebrospinal_fluid",
                  "grey_matter_in_brain_stem")

for (i in seq_along(group1_cols)) {
  data[[group1_names[i]]] <- log2_z(data[[group1_cols[i]]] * data[[head_scaling]])
}

## Group 2: volumes already normalised for head size
group2_ids   <- c("25009", "25005", "25007")
group2_cols  <- paste0("X", group2_ids, ".2.0")
group2_names <- c("grey_and_white_matter", "grey_matter", "white_matter")

for (i in seq_along(group2_cols)) {
  data[[group2_names[i]]] <- log2_z(data[[group2_cols[i]]])
}

## Group 3: left + right hemispheres summed on the raw scale, then transformed
left_ids  <- c("25878", "25880", "25882", "25884", "25886", "25888", "25890")
right_ids <- c("25879", "25881", "25883", "25885", "25887", "25889", "25891")
left_cols  <- paste0("X", left_ids,  ".2.0")
right_cols <- paste0("X", right_ids, ".2.0")
group3_names <- c("Thalamus", "Caudate_Nucleus", "Putamen", "Globus_Pallidus",
                  "Hippocampus", "Amygdala", "Nucleus_Accumbens")

for (i in seq_along(left_cols)) {
  # NA in either hemisphere propagates: a one-sided total is not comparable
  total <- (data[[left_cols[i]]] + data[[right_cols[i]]]) * data[[head_scaling]]
  data[[group3_names[i]]] <- log2_z(total)
}

data$center <- as.factor(data$center)

predictors  <- c(covariates_base, "center")
y_variables <- c(group3_names, group1_names, group2_names)

sample_counts <- colSums(!is.na(data[y_variables]))
print(range(sample_counts))

results_brain <- run_models(data, y_variables, predictors)
write.csv(results_brain, paste0(dir_out, "brain0.csv"), row.names = FALSE)


### 2. Brain MRI (cerebellar regions) -----------------------------------------

datascore <- load_score()
brain     <- read.csv(paste0(dir_raw, "comp_grey_volume_out.csv"))
center    <- read.csv(path_center)

data <- merge(datascore, brain,  by = "n_eid")
data <- merge(data,      center, by = "n_eid")
data <- merge(data,      data3,  by = "n_eid")

# Field IDs belonging to each cerebellar lobule
regions_list <- list(
  crus_i  = c("25900", "25901", "25902"),
  crus_ii = c("25903", "25904", "25905"),
  i_iv    = c("25893", "25894"),
  v       = c("25895", "25896"),
  vi      = c("25897", "25898", "25899"),
  vii     = c("25906", "25907", "25908"),
  viii    = c("25909", "25910", "25911", "25912", "25913", "25914"),
  ix      = c("25915", "25916", "25917"),
  x       = c("25918", "25919", "25920")
)

# FIX: rowSums now propagates NA instead of na.rm = TRUE. Previously a fully
# missing row summed to 0, and log2(0) = -Inf contaminated scale() and the fits.
for (region in names(regions_list)) {
  region_cols <- paste0("n_", regions_list[[region]], "_2_0")
  missing_cols <- setdiff(region_cols, colnames(data))
  if (length(missing_cols) > 0) {
    warning("Missing columns for region ", region, ": ",
            paste(missing_cols, collapse = ", "))
  }
  region_cols <- intersect(region_cols, colnames(data))
  if (length(region_cols) == 0) next
  data[[region]] <- log2_z(rowsum_complete(data[, region_cols, drop = FALSE]))
}

# FIX: center is now a factor here too (it was left numeric in the original)
data$center <- as.factor(data$center)

predictors <- c(covariates_base, "center")

# NOTE: the original listed 8 outcomes for the sample-size check but looped over
# all 9 lobules. All 9 are now modelled and FDR-corrected together. If "vii" was
# meant to be excluded, drop it from y_variables below and it will also be
# removed from the multiple-testing correction.
y_variables <- names(regions_list)

sample_counts <- colSums(!is.na(data[y_variables]))
print(range(sample_counts))

results_cereb <- run_models(data, y_variables, predictors)
write.csv(results_cereb, paste0(dir_out, "brain2.csv"), row.names = FALSE)


### 3. Vascular function ------------------------------------------------------

datascore <- load_score()
physical  <- read.csv(paste0(dir_raw, "data_update2.csv"))

data <- merge(datascore, physical, by = "n_eid")
data <- merge(data,      data3,    by = "n_eid")

# FIX: outcomes are now named explicitly. The original grepped every "_2_0"
# column, which swept in unrelated fields and fitted far more models than the
# two blood-pressure outcomes the section is about.
bp_fields <- c(diastolic_bp = "n_4079_2_0",
               systolic_bp  = "n_4080_2_0")

bp_fields <- bp_fields[bp_fields %in% colnames(data)]
if (length(bp_fields) < 2) {
  warning("Expected blood-pressure columns not all found: ",
          paste(setdiff(c("n_4079_2_0", "n_4080_2_0"), colnames(data)), collapse = ", "))
}

# FIX: the transformed values are written to NEW, explicitly named columns.
# The original created "<col>_scaled" columns but then regressed on the
# untransformed originals, because sub("_2_0", "") renamed the raw columns to
# match the names stored in scaled_vars. The scaled columns were never used.
y_variables <- names(bp_fields)
for (i in seq_along(bp_fields)) {
  data[[names(bp_fields)[i]]] <- log2_z(data[[bp_fields[i]]])
}

sample_counts <- colSums(!is.na(data[y_variables]))
print(range(sample_counts))

results_bp <- run_models(data, y_variables, covariates_base)
write.csv(results_bp, paste0(dir_out, "pressure.csv"), row.names = FALSE)


### 4. Liver biomarkers and immune cells --------------------------------------

biomarker_names <- read.csv(paste0(dir_raw, "biochemistry name.csv"))
datascore       <- load_score()
biochemistry    <- read.csv(paste0(dir_raw, "data_update.csv"))

data <- merge(datascore, data3,        by = "n_eid")
data <- merge(data,      biochemistry, by = "n_eid")

# Instance-1 biomarker measurements, taken only from the biochemistry table so
# that covariates and other merged columns cannot leak into the outcome list
biochem_cols <- grep("_1_0$", colnames(biochemistry), value = TRUE)
biochem_cols <- intersect(biochem_cols, colnames(data))

# FIX: same scaled-column bug as section 3. Transformed values go into new
# columns with a "z_" prefix, and those are what the models use.
y_variables <- paste0("z_", sub("_1_0$", "", biochem_cols))
for (i in seq_along(biochem_cols)) {
  data[[y_variables[i]]] <- log2_z(data[[biochem_cols[i]]])
}

sample_counts <- colSums(!is.na(data[y_variables]))
print(range(sample_counts))

results_bio <- run_models(data, y_variables, covariates_base)

# Map back to the original field names for the lookup table merge
results_bio$Composition <- sub("^z_", "", results_bio$Variable)
results_bio <- merge(results_bio, biomarker_names, by = "Composition", all.x = TRUE)
write.csv(results_bio, paste0(dir_out, "blood.csv"), row.names = FALSE)


### 5. Bioelectrical impedance / body composition -----------------------------

comp_names <- read.csv(paste0(dir_raw, "in_name.csv"))
datascore  <- load_score()
Dex        <- read.csv(paste0(dir_raw, "BIC_data.csv"))

data <- merge(Dex,  datascore, by = "n_eid")
data <- merge(data, data3,     by = "n_eid")

# Left/right pairs that must be combined BEFORE any transformation
paired_vars <- list(
  Leg_fat_free_mass    = c("23113", "23117"),
  Leg_fat_mass         = c("23112", "23116"),
  Leg_fat_percentage   = c("23111", "23115"),
  Leg_predicted_mass   = c("23114", "23118"),
  Arm_fat_percentage   = c("23119", "23123"),
  Arm_fat_mass         = c("23120", "23124"),
  Arm_fat_free_mass    = c("23121", "23125"),
  Arm_predicted_mass   = c("23122", "23126"),
  Impedance_arm        = c("23109", "23110"),
  Impedance_leg        = c("23107", "23108"),
  Impedance_arm_manual = c("6221",  "6222"),
  Impedance_leg_manual = c("6219",  "6220")
)

paired_ids <- unlist(paired_vars, use.names = FALSE)

# FIX: the original combined left and right AFTER differencing and
# standardising, then used rowSums(na.rm = TRUE). Adding two z-scores does not
# give a bilateral total, and na.rm = TRUE silently treated a missing limb as
# zero. Left and right are now summed on the raw scale at each visit, before
# the change score is computed.
#
# For percentage fields the sum of two limbs is not meaningful; the mean is used
# instead so the value stays on the original percentage scale.
combine_fun <- function(varname) {
  if (grepl("percentage", varname)) function(l, r) (l + r) / 2 else function(l, r) l + r
}

y_variables <- character(0)

for (new_var in names(paired_vars)) {
  ids <- paired_vars[[new_var]]
  cols_0 <- paste0("n_", ids, "_0_0")
  cols_2 <- paste0("n_", ids, "_2_0")
  
  if (!all(c(cols_0, cols_2) %in% colnames(data))) {
    warning("Missing column(s), skipping variable: ", new_var)
    next
  }
  
  f <- combine_fun(new_var)
  base_val <- f(data[[cols_0[1]]], data[[cols_0[2]]])   # NA propagates
  foll_val <- f(data[[cols_2[1]]], data[[cols_2[2]]])
  
  # FIX: change scores are standardised directly. The original shifted by the
  # sample minimum and took log2, which is data-dependent, extremely sensitive
  # to outliers, and returns Inf when a variable is entirely missing. A change
  # score is symmetric around zero and needs no log transform.
  data[[new_var]] <- as.vector(scale(foll_val - base_val))
  y_variables <- c(y_variables, new_var)
}

# Unpaired fields measured at both visits, handled the same way
all_cols <- colnames(data)
both_visit_ids <- gsub("^n_(\\d+)_[02]_0$", "\\1",
                       grep("^n_\\d+_[02]_0$", all_cols, value = TRUE))
both_visit_ids <- unique(both_visit_ids)
both_visit_ids <- setdiff(both_visit_ids, paired_ids)   # paired ones already done

for (id in both_visit_ids) {
  col_0 <- paste0("n_", id, "_0_0")
  col_2 <- paste0("n_", id, "_2_0")
  if (!all(c(col_0, col_2) %in% all_cols)) next
  
  change <- data[[col_2]] - data[[col_0]]
  if (all(is.na(change))) {
    warning("All values missing for field ", id, "; skipped.")
    next
  }
  
  new_col <- paste0("n_", id)
  data[[new_col]] <- as.vector(scale(change))
  y_variables <- c(y_variables, new_col)
}

sample_counts <- colSums(!is.na(data[y_variables]))
print(range(sample_counts))

results_imp <- run_models(data, y_variables, covariates_base)

results_imp$Composition <- results_imp$Variable
results_imp <- merge(results_imp, comp_names, by = "Composition", all.x = TRUE)

# Full results saved before filtering, so the FDR denominator stays inspectable
write.csv(results_imp, paste0(dir_out, "impedance_all.csv"), row.names = FALSE)
write.csv(results_imp %>% filter(FDR_P < 0.05),
          paste0(dir_out, "impedance.csv"), row.names = FALSE)


### 6. Cardiac function -------------------------------------------------------

datascore <- load_score()
cardiac   <- read.csv("D:/科研数据/课题数据/UKB分析/UKB数据常用汇总/MRI data/CardMRI.csv")

data <- merge(datascore, cardiac, by = "n_eid")
data <- merge(data,      data3,   by = "n_eid")

# NOTE: the original selected outcomes by column position, c(2:7, 9:84), which
# breaks if the source file's column order ever changes. Selection is now by
# name: everything except the ID and the one skipped column.
cardiac_cols <- setdiff(colnames(cardiac), "n_eid")

# Column 8 of the original file was deliberately excluded; preserved here
if (length(colnames(cardiac)) >= 8) {
  cardiac_cols <- setdiff(cardiac_cols, colnames(cardiac)[8])
}

# Keep numeric outcomes only
cardiac_cols <- cardiac_cols[sapply(data[cardiac_cols], is.numeric)]

y_variables <- paste0("z_", cardiac_cols)
for (i in seq_along(cardiac_cols)) {
  data[[y_variables[i]]] <- log2_z(data[[cardiac_cols[i]]])
}

sample_counts <- colSums(!is.na(data[y_variables]))
print(range(sample_counts))

# The original section stopped after counting samples; models are now fitted
results_cardiac <- run_models(data, y_variables, covariates_base)
results_cardiac$Variable <- sub("^z_", "", results_cardiac$Variable)
write.csv(results_cardiac, paste0(dir_out, "cardiac.csv"), row.names = FALSE)


### Optional: global FDR across all sections ----------------------------------
## Each section above corrects within itself. If a study-wide correction is
## wanted, combine and re-adjust:
# all_results <- bind_rows(
#   results_brain   %>% mutate(Domain = "Brain"),
#   results_cereb   %>% mutate(Domain = "Cerebellum"),
#   results_bp      %>% mutate(Domain = "Vascular"),
#   results_bio     %>% mutate(Domain = "Biochemistry"),
#   results_imp     %>% mutate(Domain = "BodyComposition"),
#   results_cardiac %>% mutate(Domain = "Cardiac")
# )
# all_results$FDR_P_global <- p.adjust(all_results$P_value, method = "fdr")
# write.csv(all_results, paste0(dir_out, "all_results.csv"), row.names = FALSE)