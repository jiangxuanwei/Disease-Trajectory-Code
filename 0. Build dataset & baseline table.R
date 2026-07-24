# ==============================================================================
# Table 1: Baseline characteristics of the four analytical datasets
#
# Scenario 1:
#   Training set   = BMI_2 - BMI_0
#   Validation set = BMI_1 - BMI_0, excluding participants with BMI_2
#
# Scenario 2:
#   Use BMI_2 when available; otherwise use BMI_1
#   Random 60:40 training-validation split
#
# Packages:
#   tableone + flextable + officer
#
# Only the final combined Word table is saved.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Install packages once if necessary
# ------------------------------------------------------------------------------

# install.packages(c("tableone", "flextable", "officer"))


# ------------------------------------------------------------------------------
# 2. Load packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tableone)
library(flextable)
library(officer)


# ------------------------------------------------------------------------------
# 3. File paths
# ------------------------------------------------------------------------------

file_path <- paste0(
  "D:/科研数据/课题数据/1. UKB数据课题/",
  "BMI Loss 课题/更新后结果/raw_data.csv"
)

output_dir <- paste0(
  "D:/科研数据/课题数据/1. UKB数据课题/",
  "BMI Loss 课题/更新后结果"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------------------------
# 4. Define variables
# ------------------------------------------------------------------------------

# Categorical variables
categorical_vars <- c(
  "Sex",
  "Ethnicity",
  "Edu",
  "Employed",
  "Smoke",
  "Drink"
)

# Variables included in Table 1
table1_vars <- c(
  "Age",
  "Sex",
  "Ethnicity",
  "Edu",
  "Employed",
  "Tdi",
  "Smoke",
  "Drink",
  "Mets",
  "Dietscore",
  "BMI"
)

# Variables presented as median (Q1, Q3)
# Their P values will be calculated using nonparametric tests
nonnormal_vars <- c(
  "Tdi",
  "Mets",
  "Dietscore"
)


# ------------------------------------------------------------------------------
# 5. Read and prepare the original dataset
# ------------------------------------------------------------------------------

raw_data <- read.csv(
  file_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_vars <- unique(
  c(
    "BMI_0",
    "BMI_1",
    "BMI_2",
    table1_vars
  )
)

missing_vars <- setdiff(
  required_vars,
  names(raw_data)
)

if (length(missing_vars) > 0) {
  stop(
    paste0(
      "The following variables are missing from raw_data.csv: ",
      paste(missing_vars, collapse = ", ")
    )
  )
}

# Convert categorical variables to factors before splitting.
# This helps preserve the same factor levels across the four datasets.
raw_data <- raw_data %>%
  mutate(
    across(
      all_of(categorical_vars),
      as.factor
    )
  )


# ==============================================================================
# Scenario 1: Assessment-time split
# ==============================================================================


# ------------------------------------------------------------------------------
# 6. Scenario 1 training dataset
# ------------------------------------------------------------------------------

# BMI change calculated using BMI_2 and BMI_0
scenario1_train <- raw_data %>%
  mutate(
    BMIchange = BMI_2 - BMI_0
  ) %>%
  filter(
    !is.na(BMIchange)
  ) %>%
  mutate(
    Group = factor(
      "Training",
      levels = c(
        "Training",
        "Validation"
      )
    )
  )


# ------------------------------------------------------------------------------
# 7. Scenario 1 validation dataset
# ------------------------------------------------------------------------------

# BMI change calculated using BMI_1 and BMI_0
# Participants with BMI_2 are excluded to avoid overlap
scenario1_validation <- raw_data %>%
  mutate(
    BMIchange = BMI_1 - BMI_0
  ) %>%
  filter(
    !is.na(BMIchange),
    is.na(BMI_2)
  ) %>%
  mutate(
    Group = factor(
      "Validation",
      levels = c(
        "Training",
        "Validation"
      )
    )
  )


# Combine the Scenario 1 datasets
scenario1_data <- bind_rows(
  scenario1_train,
  scenario1_validation
)


# ==============================================================================
# Scenario 2: Random 60:40 split
# ==============================================================================


# ------------------------------------------------------------------------------
# 8. Construct the overall Scenario 2 dataset
# ------------------------------------------------------------------------------

# Use BMI_2 when available; otherwise use BMI_1
scenario2_all <- raw_data %>%
  mutate(
    BMIchange = case_when(
      !is.na(BMI_0) & !is.na(BMI_2) ~
        BMI_2 - BMI_0,
      
      !is.na(BMI_0) &
        is.na(BMI_2) &
        !is.na(BMI_1) ~
        BMI_1 - BMI_0,
      
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(BMIchange)
  )


# ------------------------------------------------------------------------------
# 9. Apply the same 60:40 random split as the main analysis
# ------------------------------------------------------------------------------

set.seed(123)

n_total <- nrow(scenario2_all)

train_index <- sample(
  seq_len(n_total),
  size = round(0.6 * n_total)
)


# Scenario 2 training dataset
scenario2_train <- scenario2_all[
  train_index,
  ,
  drop = FALSE
] %>%
  mutate(
    Group = factor(
      "Training",
      levels = c(
        "Training",
        "Validation"
      )
    )
  )


# Scenario 2 validation dataset
scenario2_validation <- scenario2_all[
  -train_index,
  ,
  drop = FALSE
] %>%
  mutate(
    Group = factor(
      "Validation",
      levels = c(
        "Training",
        "Validation"
      )
    )
  )


# Combine the Scenario 2 datasets
scenario2_data <- bind_rows(
  scenario2_train,
  scenario2_validation
)


# ------------------------------------------------------------------------------
# 10. Check sample sizes
# ------------------------------------------------------------------------------

cat(
  "\nScenario 1 training N = ",
  nrow(scenario1_train),
  "\n",
  sep = ""
)

cat(
  "Scenario 1 validation N = ",
  nrow(scenario1_validation),
  "\n",
  sep = ""
)

cat(
  "Scenario 2 training N = ",
  nrow(scenario2_train),
  "\n",
  sep = ""
)

cat(
  "Scenario 2 validation N = ",
  nrow(scenario2_validation),
  "\n\n",
  sep = ""
)


# ==============================================================================
# Generate Table 1 using tableone
# ==============================================================================


# ------------------------------------------------------------------------------
# 11. Function to create and extract each Table 1
# ------------------------------------------------------------------------------

make_table1 <- function(data) {
  
  # Create TableOne object
  table_object <- CreateTableOne(
    vars = table1_vars,
    strata = "Group",
    data = data,
    factorVars = categorical_vars,
    includeNA = FALSE,
    test = TRUE,
    smd = FALSE,
    addOverall = FALSE
  )
  
  # Convert TableOne output to a regular matrix
  #
  # Continuous variables not included in nonnormal_vars:
  #   mean (SD), ordinary parametric comparison
  #
  # Variables included in nonnormal_vars:
  #   median (Q1, Q3), nonparametric comparison
  #
  # Categorical variables:
  #   n (%), chi-square test
  table_matrix <- print(
    table_object,
    showAllLevels = TRUE,
    nonnormal = nonnormal_vars,
    catDigits = 1,
    contDigits = 2,
    pDigits = 3,
    quote = FALSE,
    missing = FALSE,
    explain = TRUE,
    test = TRUE,
    smd = FALSE,
    noSpaces = TRUE,
    printToggle = FALSE,
    formatOptions = list(
      scientific = FALSE,
      big.mark = ","
    )
  )
  
  # Convert to data frame
  table_df <- as.data.frame(
    table_matrix,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Move row names into the first column
  table_df <- data.frame(
    Characteristic = rownames(table_matrix),
    table_df,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  rownames(table_df) <- NULL
  
  # Remove spaces from column names
  names(table_df) <- trimws(
    names(table_df)
  )
  
  # Remove the test-name column, if generated
  test_column <- which(
    tolower(names(table_df)) == "test"
  )
  
  if (length(test_column) > 0) {
    table_df <- table_df[
      ,
      -test_column,
      drop = FALSE
    ]
  }
  
  # Identify the level column
  level_column <- names(table_df)[
    tolower(names(table_df)) == "level"
  ]
  
  if (length(level_column) == 0) {
    table_df$level <- ""
    level_column <- "level"
  }
  
  # Identify the P-value column
  p_column <- names(table_df)[
    tolower(names(table_df)) == "p"
  ]
  
  if (length(p_column) == 0) {
    stop(
      "The P-value column could not be identified in the tableone output."
    )
  }
  
  # Check group columns
  if (
    !"Training" %in% names(table_df) ||
    !"Validation" %in% names(table_df)
  ) {
    stop(
      paste0(
        "Training or Validation column could not be identified. ",
        "Current column names are: ",
        paste(names(table_df), collapse = ", ")
      )
    )
  }
  
  # Keep and standardize the required columns
  output_df <- data.frame(
    Characteristic = table_df$Characteristic,
    Level = table_df[[level_column[1]]],
    Training = table_df[["Training"]],
    Validation = table_df[["Validation"]],
    P_value = table_df[[p_column[1]]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  output_df$Characteristic <- trimws(
    output_df$Characteristic
  )
  
  output_df$Level <- ifelse(
    is.na(output_df$Level),
    "",
    trimws(output_df$Level)
  )
  
  return(output_df)
}


# ------------------------------------------------------------------------------
# 12. Generate Scenario 1 and Scenario 2 tables
# ------------------------------------------------------------------------------

table1_scenario1 <- make_table1(
  scenario1_data
)

table1_scenario2 <- make_table1(
  scenario2_data
)


# ------------------------------------------------------------------------------
# 13. Apply publication-style variable labels
# ------------------------------------------------------------------------------

label_map <- c(
  "Age" = "Age, years",
  "Sex" = "Sex",
  "Ethnicity" = "Ethnicity",
  "Edu" = "Education",
  "Employed" = "Employment status",
  "Tdi" = "Townsend deprivation index",
  "Smoke" = "Smoking status",
  "Drink" = "Drinking status",
  "Mets" = "Physical activity",
  "Dietscore" = "Diet score",
  "BMI" = "Baseline BMI, kg/m²"
)


replace_variable_labels <- function(x) {
  
  for (old_name in names(label_map)) {
    
    pattern <- paste0(
      "^",
      old_name,
      "(?=\\s|\\(|$)"
    )
    
    x <- sub(
      pattern,
      label_map[[old_name]],
      x,
      perl = TRUE
    )
  }
  
  return(x)
}


table1_scenario1$Characteristic <- replace_variable_labels(
  table1_scenario1$Characteristic
)

table1_scenario2$Characteristic <- replace_variable_labels(
  table1_scenario2$Characteristic
)


# ------------------------------------------------------------------------------
# 14. Create row identifiers to align both scenario tables
# ------------------------------------------------------------------------------

create_row_key <- function(table_df) {
  
  characteristic_filled <- table_df$Characteristic
  
  for (i in seq_along(characteristic_filled)) {
    
    if (
      i > 1 &&
      (
        is.na(characteristic_filled[i]) ||
        characteristic_filled[i] == ""
      )
    ) {
      characteristic_filled[i] <- characteristic_filled[i - 1]
    }
  }
  
  paste(
    characteristic_filled,
    table_df$Level,
    sep = "|||"
  )
}


scenario1_key <- create_row_key(
  table1_scenario1
)

scenario2_key <- create_row_key(
  table1_scenario2
)

scenario2_order <- match(
  scenario1_key,
  scenario2_key
)

if (any(is.na(scenario2_order))) {
  stop(
    paste0(
      "The rows of the two scenario tables could not be fully aligned. ",
      "Please check whether factor levels are consistent across datasets."
    )
  )
}


# ------------------------------------------------------------------------------
# 15. Combine the two scenarios
# ------------------------------------------------------------------------------

table1_combined <- data.frame(
  Characteristic = table1_scenario1$Characteristic,
  Level = table1_scenario1$Level,
  
  S1_Training = table1_scenario1$Training,
  S1_Validation = table1_scenario1$Validation,
  S1_P = table1_scenario1$P_value,
  
  S2_Training =
    table1_scenario2$Training[scenario2_order],
  
  S2_Validation =
    table1_scenario2$Validation[scenario2_order],
  
  S2_P =
    table1_scenario2$P_value[scenario2_order],
  
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Display the combined data frame in R
print(
  table1_combined,
  row.names = FALSE
)


# ==============================================================================
# Format and save the final combined Word table
# ==============================================================================


# ------------------------------------------------------------------------------
# 16. Create flextable
# ------------------------------------------------------------------------------

table1_flex <- flextable(
  table1_combined
)


# Create two-level column headers
header_mapping <- data.frame(
  col_keys = names(table1_combined),
  
  Header1 = c(
    "Characteristic",
    "Level",
    rep(
      "Scenario 1: assessment-time split",
      3
    ),
    rep(
      "Scenario 2: 60:40 random split",
      3
    )
  ),
  
  Header2 = c(
    "Characteristic",
    "Level",
    "Training",
    "Validation",
    "P value",
    "Training",
    "Validation",
    "P value"
  ),
  
  stringsAsFactors = FALSE
)


table1_flex <- table1_flex %>%
  
  set_header_df(
    mapping = header_mapping,
    key = "col_keys"
  ) %>%
  
  # Merge the scenario headings horizontally
  merge_h(
    part = "header"
  ) %>%
  
  # Merge Characteristic and Level headings vertically
  merge_v(
    j = c(
      "Characteristic",
      "Level"
    ),
    part = "header"
  ) %>%
  
  # Format font
  font(
    fontname = "Arial",
    part = "all"
  ) %>%
  
  fontsize(
    size = 8,
    part = "all"
  ) %>%
  
  # Bold header
  bold(
    part = "header"
  ) %>%
  
  # Header alignment
  align(
    align = "center",
    part = "header"
  ) %>%
  
  # Left-align characteristic and level columns
  align(
    j = c(
      "Characteristic",
      "Level"
    ),
    align = "left",
    part = "body"
  ) %>%
  
  # Center-align numerical columns
  align(
    j = c(
      "S1_Training",
      "S1_Validation",
      "S1_P",
      "S2_Training",
      "S2_Validation",
      "S2_P"
    ),
    align = "center",
    part = "body"
  ) %>%
  
  valign(
    valign = "center",
    part = "all"
  ) %>%
  
  autofit()


# Set selected column widths
table1_flex <- table1_flex %>%
  width(
    j = "Characteristic",
    width = 2.3
  ) %>%
  width(
    j = "Level",
    width = 1.0
  )


# ------------------------------------------------------------------------------
# 17. Landscape Word page settings
# ------------------------------------------------------------------------------

landscape_section <- prop_section(
  page_size = page_size(
    orient = "landscape"
  ),
  page_margins = page_mar(
    top = 0.5,
    bottom = 0.5,
    left = 0.5,
    right = 0.5
  ),
  type = "continuous"
)


# ------------------------------------------------------------------------------
# 18. Save only the final combined table
# ------------------------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "Table1_four_datasets_combined.docx"
)

save_as_docx(
  `Table 1. Baseline characteristics of the analytical datasets` =
    table1_flex,
  
  path = output_file,
  
  pr_section = landscape_section
)


cat(
  "\nThe final combined Table 1 has been saved to:\n",
  output_file,
  "\n",
  sep = ""
)