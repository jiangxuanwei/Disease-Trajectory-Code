# ============================================================
# Subgroup Definition
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

### 1. Tumor stage------------------------------------------------------------------------

# Load the necessary datasets
track <- read.csv("D:/CMD_CA.csv")
disease <- read.csv("D:/track.csv")
covariate <- read.csv("D:/covariates_impute.csv")

# Merge the datasets by 'n_eid'
data <- merge(track, disease, by = "n_eid")
data <- merge(data, covariate, by = "n_eid")

### Tumor stage1: Filter data for tumor stages 1 and 2
data1 <- data %>% filter(event_CA == 0 | event_CA == 1 & tumorstage %in% c(1, 2))

### Tumor stage3: Filter data for tumor stage 3
data2 <- data %>% filter(event_CA == 0 | event_CA == 1 & tumorstage == 3)

write.csv(data1, "D:/CMD-CA-CA-stage1.csv")
write.csv(data2, "D:/CMD-CA-CA-stage2.csv")

### 2. CMD stage------------------------------------------------------------------------

# Load the datasets
track <- read.csv("D:/CMD_CA.csv")
disease <- read.csv("D:/track.csv")
covariate <- read.csv("D:/covariates_impute.csv")

# Merge the datasets
data <- merge(track, disease, by = "n_eid")
data <- merge(data, covariate, by = "n_eid")

# Create new variables to indicate whether certain health conditions are above certain thresholds
data$event_CD = 0
data$event_CD[data$event_dia == 1 | data$event_CVD == 1] = 1
data$time_CD = pmin(data$time_dia, data$time_CVD, na.rm = TRUE)

data <- data %>%
  mutate(BMI_status = ifelse(BMI < 30, 0, 1),
         BP_status = ifelse(sbp < 130, 0, 1),
         GLU_status = ifelse(glucose < 7, 0, 1),
         LDL_status = ifelse(ldl < 2.6, 0, 1),
         all_status = BMI_status + BP_status + GLU_status + LDL_status)

### CMD stage1: Filter data for CMD stage 1
data1 <- data %>% filter(event_CD == 0 | event_CD == 1 & all_status < 3)

### CMD stage2: Filter data for CMD stage 2
data2 <- data %>% filter(event_CD == 0 | event_CD == 1 & all_status >= 3)

write.csv(data1, "D:/CMD-CA-CMD-stage1.csv")
write.csv(data2, "D:/CMD-CA-CMD-stage2.csv")

### 3. Tumor and CMD stage------------------------------------------------------------------------

# Load the datasets
track <- read.csv("D:/CMD_CA.csv")
disease <- read.csv("D:/track.csv")
covariate <- read.csv("D:/covariates_impute.csv")

# Merge the datasets
data <- merge(track, disease, by = "n_eid")
data <- merge(data, covariate, by = "n_eid")

# Create new variables to indicate certain health conditions
data$event_CD = 0
data$event_CD[data$event_dia == 1 | data$event_CVD == 1] = 1
data$time_CD = pmin(data$time_dia, data$time_CVD, na.rm = TRUE)

data <- data %>%
  mutate(BMI_status = ifelse(BMI < 30, 0, 1),
         BP_status = ifelse(sbp < 130, 0, 1),
         GLU_status = ifelse(glucose < 7, 0, 1),
         LDL_status = ifelse(ldl < 2.6, 0, 1),
         all_status = BMI_status + BP_status + GLU_status + LDL_status)

# Filter the data for CMD stage1 and CA stage1, and CMD stage3 and CA stage3
data1 <- data %>% filter(event_CD == 0 | event_CD == 1 & all_status < 3) %>% filter(event_CA == 0 | event_CA == 1 & tumorstage %in% c(1, 2))
data2 <- data %>% filter(event_CD == 0 | event_CD == 1 & all_status >= 3) %>% filter(event_CA == 0 | event_CA == 1 & tumorstage == 3)

write.csv(data1, "D:/CMD-CA-both-stage1.csv")
write.csv(data2, "D:/CMD-CA-both-stage3.csv")

### Cancer categories------------------------------------------------------------------------

# Load the datasets
track <- read.csv("D:/CMD_CA.csv")
disease <- read.csv("D:/track.csv")
covariate <- read.csv("D:/covariates_impute.csv")

# Merge the datasets
data <- merge(track, disease, by = "n_eid")
data <- merge(data, covariate, by = "n_eid")

### High survival rate cancer: Filter data for cancers with high survival rates
data1 <- data %>% filter(event_CA == 0 | event_SkinCA == 1 | event_ProsCA == 1 | event_BreastCA == 1)

### Low survival rate cancer: Filter data for cancers with low survival rates
data2 <- data %>% filter(event_CA == 0 | event_LiverCA == 1 | event_LungCA == 1 |
                           event_EsophagealCA == 1 | event_PancreaticCA == 1)

write.csv(data1, "D:/CMD-CA-CA-cat1.csv")
write.csv(data2, "D:/CMD-CA-CA-cat2.csv")
