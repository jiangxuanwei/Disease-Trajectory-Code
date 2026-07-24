library(dplyr)
library(survival)
library(glmnet)
library(ggplot2)
library(table1)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

raw_data_file <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/更新后结果/raw_data.csv"
lasso_coef_file <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/更新后结果/LASSO_coefficients_decrease.csv"
bmi_score_file <- "D:/科研数据/课题数据/UKB分析/BMI Loss 课题/更新后结果/BMIscore_predict.csv"

categorical_vars <- c("Sex", "Ethnicity", "Edu", "Employed", "Smoke", "Drink", "Mets")

# ------------------------------------------------------------------------------
# 1. Build BMI loss score using the training set
# ------------------------------------------------------------------------------

set.seed(123)
raw_data <- read.csv(raw_data_file)
raw_data$BMIchange <- raw_data$BMI_2 - raw_data$BMI_0

data <- raw_data %>% filter(!is.na(BMIchange))

# Candidate proteins used for LASSO model construction.
proteins <- c("ache", "adamts8", "angptl7", "cd93", "comp", "dlk1", "fap", "gcg", "gfra3", "igsf3", "ism1", "itgb6", "prtg", "ret", "thbs4")

x_data <- data %>% select(all_of(proteins))
y_data <- data$BMIchange
x_matrix <- as.matrix(x_data)
y_vector <- as.vector(y_data)

# Fit LASSO model using cross-validation.
lasso_model <- cv.glmnet(x_matrix, y_vector, alpha = 1)
lasso_model$lambda.min

# Extract LASSO coefficients at lambda.min.
lasso_coefficients <- coef(lasso_model, s = "lambda.min")
lasso_coefficients_df <- as.data.frame(as.matrix(lasso_coefficients))
write.csv(lasso_coefficients_df, lasso_coef_file)

# Calculate BMI loss prediction score.
coef_de <- read.csv(lasso_coef_file)
variables <- coef_de$Variable
coefficients <- coef_de$Coefficient
matched_vars <- intersect(variables, colnames(raw_data))

BMI_predict_matrix <- sweep(raw_data[, matched_vars, drop = FALSE], 2, coefficients, `*`)
raw_data$BMI_predict_decrease <- rowSums(BMI_predict_matrix)

datascore <- raw_data[, c("n_eid", "BMI_predict_decrease")]
write.csv(datascore, bmi_score_file)


# ------------------------------------------------------------------------------
# 2. Analysis stratified by follow-up assessment
# ------------------------------------------------------------------------------

# Correlation analysis: training set.
raw_data <- read.csv(raw_data_file)
datascore <- read.csv(bmi_score_file)
datascore$BMI_predict_decrease <- datascore$BMI_predict_decrease * (-1)

data <- merge(raw_data, datascore, by = "n_eid")
data$BMIchange <- data$BMI_2 - data$BMI_0
data <- data %>% filter(!is.na(BMIchange))

ggplot(data, aes(x = BMI_predict_decrease, y = BMIchange)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_y_continuous(limits = c(-20, 20)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  labs(x = "Predicted BMI Change", y = "Observed BMI Change") +
  theme_minimal()

# Correlation analysis: validation set.
data <- merge(raw_data, datascore, by = "n_eid")
data$BMIchange <- data$BMI_1 - data$BMI_0
data <- data %>% filter(!is.na(BMIchange)) %>% filter(is.na(BMI_2))

ggplot(data, aes(x = BMI_predict_decrease, y = BMIchange)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_y_continuous(limits = c(-20, 20)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  labs(x = "Predicted BMI Change", y = "Observed BMI Change") +
  theme_minimal()

# Linear association: training set.
data <- merge(raw_data, datascore, by = "n_eid")
data$BMIchange <- data$BMI_2 - data$BMI_0
data <- data %>% filter(!is.na(BMIchange)) %>% mutate(across(all_of(categorical_vars), as.factor))
data$score <- ntile(data$BMI_predict_decrease, 3)

table1(~ BMIchange | score, data)

a <- lm(BMIchange ~ BMI_predict_decrease + BMI_0 + Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore, data = data)
summary_a_continuous <- summary(a)

a <- lm(BMIchange ~ as.factor(score) + BMI_0 + Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore, data = data)
summary_a_tertile <- summary(a)

# Linear association: validation set.
data <- merge(raw_data, datascore, by = "n_eid")
data$BMIchange <- data$BMI_1 - data$BMI_0
data <- data %>% filter(!is.na(BMIchange)) %>% filter(is.na(BMI_2)) %>% mutate(across(all_of(categorical_vars), as.factor))
data$score <- ntile(data$BMI_predict_decrease, 3)

table1(~ BMIchange | score, data)

a <- lm(BMIchange ~ BMI_predict_decrease + BMI_0 + Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore, data = data)
summary_a_continuous <- summary(a)

a <- lm(BMIchange ~ as.factor(score) + BMI_0 + Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore, data = data)
summary_a_tertile <- summary(a)


# ------------------------------------------------------------------------------
# 3. Analysis using a 60:40 train-test split
# ------------------------------------------------------------------------------

raw_data <- read.csv(raw_data_file)
datascore <- read.csv(bmi_score_file)
datascore$BMI_predict_decrease <- datascore$BMI_predict_decrease * (-1)

data <- merge(raw_data, datascore, by = "n_eid")
data <- data %>% mutate(BMIchange = if_else(is.na(BMI_2), BMI_1 - BMI_0, BMI_2 - BMI_0))
data <- data %>% filter(!is.na(BMIchange)) %>% mutate(across(all_of(categorical_vars), as.factor))

set.seed(123)
n <- nrow(data)
train_index <- sample(seq_len(n), size = round(0.6 * n))
train_data <- data[train_index, ]
test_data <- data[-train_index, ]

# Correlation analysis: training set.
ggplot(train_data, aes(x = BMI_predict_decrease, y = BMIchange)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_y_continuous(limits = c(-20, 20)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  labs(x = "Predicted BMI Change", y = "Observed BMI Change") +
  theme_minimal()

# Correlation analysis: validation set.
ggplot(test_data, aes(x = BMI_predict_decrease, y = BMIchange)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_y_continuous(limits = c(-20, 20)) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  labs(x = "Predicted BMI Change", y = "Observed BMI Change") +
  theme_minimal()

# Linear association: training set.
train_data$score <- ntile(train_data$BMI_predict_decrease, 3)

table1(~ BMIchange | score, train_data)

a <- lm(BMIchange ~ BMI_predict_decrease + BMI_0 + Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore, data = train_data)
summary_a_continuous <- summary(a)

a <- lm(BMIchange ~ as.factor(score) + BMI_0 + Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore, data = train_data)
summary_a_tertile <- summary(a)

# Linear association: validation set.
test_data$score <- ntile(test_data$BMI_predict_decrease, 3)

table1(~ BMIchange | score, test_data)

a <- lm(BMIchange ~ BMI_predict_decrease + BMI_0 + Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore, data = test_data)
summary_a_continuous <- summary(a)

a <- lm(BMIchange ~ as.factor(score) + BMI_0 + Age + Sex + Ethnicity + Edu + Employed + Tdi + Smoke + Drink + Mets + Dietscore, data = test_data)
summary_a_tertile <- summary(a)