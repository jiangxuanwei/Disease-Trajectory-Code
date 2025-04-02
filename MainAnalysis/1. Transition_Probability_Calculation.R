# ============================================================
# UK Biobank Disease Transition Probability Analysis
# Author: JXW
# Last Modified: Feb 2025
# 
# This script analyzes transition probabilities between:
# 1. Baseline to first disease (CMD/CA/Death)
# 2. CMD to multimorbidity or mortality
# 3. CA to multimorbidity or mortality 
# 4. CMD-CA multimorbidity prognosis
# 5. CA-CMD multimorbidity prognosis
# ============================================================

### 1. Baseline to 1st disease-------------------------------
# Load the dplyr library
library(dplyr)
# Read the track and covariate data from CSV files
track<-read.csv("CMD_CA.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
# Select relevant columns from the dataset for each transition pattern
datatrans1<-data[,c(1,2,3,4,30:47)]
datatrans2<-data[,c(1,5,6,7,30:47)]
datatrans3<-data[,c(1,8,9,10,30:47)]
# Assign transition patterns to each dataset
datatrans1$transpattern<-1
datatrans2$transpattern<-2
datatrans3$transpattern<-3
# Rename the columns for each dataset
names(datatrans1)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
names(datatrans2)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
names(datatrans3)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
dataall1<-rbind(datatrans1,datatrans2,datatrans3)
# Calculate the time to the next event
dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
library("mstate")
# Define the transition matrix for multi-state analysis (Baseline -> CMD, CA, Death)
tmat <- transMat(x = list(c(2,3,4), c(), c(),c()),
                 names = c( "Baseline","CMD","CA","Death"))
# Define the covariates for the model
covs_P <- c("Age", "Sex", "Ethnicity", "Edu", 
            "Employed", "Tdi", "Smoke", "Drink", 
            "Mets", "Dietscore", "Sleeptime", "BMI",
            "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
            "CMD_treat", "cancer_screen")
# Expand the covariates for the dataset
msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
# Get the column names of the expanded dataset
variables <- colnames(msebmt_allname[,24:77])
# Construct the formula
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
# Calculate the mean values of the covariates
mean_values <- colMeans(data[, c("Age", "Sex", "Ethnicity", "Edu", 
                                 "Employed", "Tdi", "Smoke", "Drink", 
                                 "Mets", "Dietscore", "Sleeptime", "BMI",
                                 "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
                                 "CMD_treat", "cancer_screen")], na.rm = TRUE)
# Create a new dataset using the mean values of covariates
newd <- data.frame(
  Age = rep(mean_values["Age"], 3),
  Sex = rep(mean_values["Sex"], 3),
  Ethnicity = rep(mean_values["Ethnicity"], 3),
  Edu = rep(mean_values["Edu"], 3),
  Employed = rep(mean_values["Employed"], 3),
  Tdi = rep(mean_values["Tdi"], 3),
  Smoke = rep(mean_values["Smoke"], 3),
  Drink = rep(mean_values["Drink"], 3),
  Mets = rep(mean_values["Mets"], 3),
  Dietscore = rep(mean_values["Dietscore"], 3),
  Sleeptime = rep(mean_values["Sleeptime"], 3),
  BMI = rep(mean_values["BMI"], 3),
  T2dhis = rep(mean_values["T2dhis"], 3),
  Cvdhis = rep(mean_values["Cvdhis"], 3),
  hiscancer = rep(mean_values["hiscancer"], 3),
  cancer_treat = rep(mean_values["cancer_treat"], 3),
  CMD_treat = rep(mean_values["CMD_treat"], 3),
  cancer_screen = rep(mean_values["cancer_screen"], 3),
  trans = 1:3
)
# Assign the transition matrix to the new dataset
attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
covs_P1 <- c("Age", "Sex", "Ethnicity", "Edu", 
             "Employed", "Tdi", "Smoke", "Drink", 
             "Mets", "Dietscore", "Sleeptime", "BMI",
             "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
             "CMD_treat", "cancer_screen")
# Expand the covariates for the new dataset
newd <- expand.covs(newd, covs_P1, longnames = FALSE)
# Define the strata for the transitions
newd$strata<-1:3
# Compute the multi-state probabilities 
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
# Plot the survival curves
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)
# Estimate the mean values of probability and SE at 5, 10, 15 years
prob1_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob1_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob1_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
prob2_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper3),3),")")
prob2_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper3),3),")")
prob2_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper3),3),")")
prob3_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate4),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower4),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper4),3),")")
prob3_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate4),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower4),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper4),3),")")
prob3_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate4),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower4),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper4),3),")")
se1_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se1_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se1_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)
se2_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se3),3)
se2_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se3),3)
se2_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se3),3)
se3_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se4),3)
se3_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se4),3)
se3_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se4),3)

### 2. CMD to multimorbidity or mortality-------------------------------
track<-read.csv("CMD_CA.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
# Select individuals with one morbidity and calculate mean values of the covariates
values<- data %>% filter (trans1==1|trans2==1)
mean_values <- colMeans(values[, c("Age", "Sex", "Ethnicity", "Edu", 
                                   "Employed", "Tdi", "Smoke", "Drink", 
                                   "Mets", "Dietscore", "Sleeptime", "BMI",
                                   "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
                                   "CMD_treat", "cancer_screen")], na.rm = TRUE)
# CMD first
data<-data %>% filter (trans1==1)
# Select relevant columns from the dataset for each transition pattern
datatrans4<-data[,c(1,11,12,13,30:47)]
datatrans5<-data[,c(1,14,15,16,30:47)]
# Assign transition patterns to each dataset
datatrans4$transpattern<-1
datatrans5$transpattern<-2
# Rename the columns for each dataset
names(datatrans4)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
names(datatrans5)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
dataall1<-rbind(datatrans4,datatrans5)
# Calculate the time to the next event
dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
library("mstate")
# Define the transition matrix for multi-state analysis (CMD -> CMD-CA, Death)
tmat <- transMat(x = list(c(2,3), c(), c()),
                 names = c( "CMD","CMD-CA","Death"))
# Define the covariates for the model
covs_P <- c("Age", "Sex", "Ethnicity", "Edu", 
            "Employed", "Tdi", "Smoke", "Drink", 
            "Mets", "Dietscore", "Sleeptime", "BMI",
            "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
            "CMD_treat", "cancer_screen")
# Expand the covariates for the dataset
msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
# Get the column names of the expanded dataset
variables <- colnames(msebmt_allname[,24:59])
# Construct the formula
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
# Create a new dataset using the mean values of covariates
newd <- data.frame(
  Age = rep(mean_values["Age"], 2),
  Sex = rep(mean_values["Sex"], 2),
  Ethnicity = rep(mean_values["Ethnicity"], 2),
  Edu = rep(mean_values["Edu"], 2),
  Employed = rep(mean_values["Employed"], 2),
  Tdi = rep(mean_values["Tdi"], 2),
  Smoke = rep(mean_values["Smoke"], 2),
  Drink = rep(mean_values["Drink"], 2),
  Mets = rep(mean_values["Mets"], 2),
  Dietscore = rep(mean_values["Dietscore"], 2),
  Sleeptime = rep(mean_values["Sleeptime"], 2),
  BMI = rep(mean_values["BMI"], 2),
  T2dhis = rep(mean_values["T2dhis"], 2),
  Cvdhis = rep(mean_values["Cvdhis"], 2),
  hiscancer = rep(mean_values["hiscancer"], 2),
  cancer_treat = rep(mean_values["cancer_treat"], 2),
  CMD_treat = rep(mean_values["CMD_treat"], 2),
  cancer_screen = rep(mean_values["cancer_screen"], 2),
  trans = 1:2
)
# Assign the transition matrix to the new dataset
attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
covs_P1 <- c("Age", "Sex", "Ethnicity", "Edu", 
             "Employed", "Tdi", "Smoke", "Drink", 
             "Mets", "Dietscore", "Sleeptime", "BMI",
             "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
             "CMD_treat", "cancer_screen")
# Expand the covariates for the new dataset
newd <- expand.covs(newd, covs_P1, longnames = FALSE)
# Define the strata for the transitions
newd$strata<-1:2
# Compute the multi-state probabilitie
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
# Plot the survival curves
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)
# Estimate the mean values of probability and standard errors at 5, 10, 15 years
prob4_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob4_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob4_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
prob5_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper3),3),")")
prob5_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper3),3),")")
prob5_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper3),3),")")
se4_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se4_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se4_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)
se5_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se3),3)
se5_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se3),3)
se5_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se3),3)

### 3. CA to multimorbidity or mortality-------------------------------
track<-read.csv("CMD_CA.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
# CA first
data<-data %>% filter (trans2==1)
# Select relevant columns from the dataset for each transition pattern
datatrans6<-data[,c(1,17,18,19,30:47)]
datatrans7<-data[,c(1,20,21,22,30:47)]
# Assign transition patterns to each dataset
datatrans6$transpattern<-1
datatrans7$transpattern<-2
# Rename the columns for each dataset
names(datatrans6)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
names(datatrans7)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
# Calculate the time to the next event
dataall1<-rbind(datatrans6,datatrans7)
dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
library("mstate")
# Define the transition matrix for multi-state analysis (CA -> CA-CMD, Death)
tmat <- transMat(x = list(c(2,3), c(), c()),
                 names = c( "CA","CA-CMD","Death"))
# Define the covariates for the model
covs_P <- c("Age", "Sex", "Ethnicity", "Edu", 
            "Employed", "Tdi", "Smoke", "Drink", 
            "Mets", "Dietscore", "Sleeptime", "BMI",
            "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
            "CMD_treat", "cancer_screen")
# Expand the covariates for the dataset
msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
# Get the column names of the expanded dataset
variables <- colnames(msebmt_allname[,24:59])
# Construct the formula
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
# Create a new dataset using the mean values of covariates
newd <- data.frame(
  Age = rep(mean_values["Age"], 2),
  Sex = rep(mean_values["Sex"], 2),
  Ethnicity = rep(mean_values["Ethnicity"], 2),
  Edu = rep(mean_values["Edu"], 2),
  Employed = rep(mean_values["Employed"], 2),
  Tdi = rep(mean_values["Tdi"], 2),
  Smoke = rep(mean_values["Smoke"], 2),
  Drink = rep(mean_values["Drink"], 2),
  Mets = rep(mean_values["Mets"], 2),
  Dietscore = rep(mean_values["Dietscore"], 2),
  Sleeptime = rep(mean_values["Sleeptime"], 2),
  BMI = rep(mean_values["BMI"], 2),
  T2dhis = rep(mean_values["T2dhis"], 2),
  Cvdhis = rep(mean_values["Cvdhis"], 2),
  hiscancer = rep(mean_values["hiscancer"], 2),
  cancer_treat = rep(mean_values["cancer_treat"], 2),
  CMD_treat = rep(mean_values["CMD_treat"], 2),
  cancer_screen = rep(mean_values["cancer_screen"], 2),
  trans = 1:2
)

# Assign the transition matrix to the new dataset
attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
covs_P1 <- c("Age", "Sex", "Ethnicity", "Edu", 
             "Employed", "Tdi", "Smoke", "Drink", 
             "Mets", "Dietscore", "Sleeptime", "BMI",
             "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
             "CMD_treat", "cancer_screen")
# Expand the covariates for the new dataset
newd <- expand.covs(newd, covs_P1, longnames = FALSE)
# Define the strata for the transitions
newd$strata<-1:2
# Compute the multi-state probabilities
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
# Plot the survival curves
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)
# Estimate the mean values of probability and standard errors at 5, 10, 15 years
prob6_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob6_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob6_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
prob7_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper3),3),")")
prob7_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper3),3),")")
prob7_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate3),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower3),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper3),3),")")
se6_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se6_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se6_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)
se7_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se3),3)
se7_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se3),3)
se7_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se3),3)

### 4. CMD-CA multimorbidity prognosis-------------------------------
track<-read.csv("CMD_CA.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
# Select individuals with multimorbidity and calculate mean values of the covariates
values<- data %>% filter (trans4==1|trans6==1)
mean_values <- colMeans(values[, c("Age", "Sex", "Ethnicity", "Edu", 
                                   "Employed", "Tdi", "Smoke", "Drink", 
                                   "Mets", "Dietscore", "Sleeptime", "BMI",
                                   "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
                                   "CMD_treat", "cancer_screen")], na.rm = TRUE)
# CMD-CA individuals
data <-data %>% filter (trans4==1)
# Define no-death outcome
data$trans10[data$trans8==1]=0
data$trans10[data$trans8==0]=1
# Select relevant columns from the dataset for each transition pattern
datatrans8<-data[,c(1,23,24,25,30:47)]
datatrans10<-data[,c(1,58,24,25,30:47)]
# Assign transition patterns to each dataset
datatrans8$transpattern<-1
datatrans10$transpattern<-2
# Rename the columns for each dataset
names(datatrans8)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
names(datatrans10)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
dataall1<-rbind(datatrans8,datatrans10)
# Calculate the time to the next event
dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
library("mstate")
# Define the transition matrix for multi-state analysis (CMD-CA -> Death, non-Death)
tmat <- transMat(x = list(c(2,3), c(), c()),
                 names = c( "CMD-CA","Death","Non-Death"))
# Define the covariates for the model
covs_P <- c("Age", "Sex", "Ethnicity", "Edu", 
            "Employed", "Tdi", "Smoke", "Drink", 
            "Mets", "Dietscore", "Sleeptime", "BMI",
            "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
            "CMD_treat", "cancer_screen")
# Expand the covariates for the dataset
msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
# Get the column names of the expanded dataset
variables <- colnames(msebmt_allname[,24:59])
# Construct the formula
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
# Create a new dataset using the mean values of covariates
newd <- data.frame(
  Age = rep(mean_values["Age"], 2),
  Sex = rep(mean_values["Sex"], 2),
  Ethnicity = rep(mean_values["Ethnicity"], 2),
  Edu = rep(mean_values["Edu"], 2),
  Employed = rep(mean_values["Employed"], 2),
  Tdi = rep(mean_values["Tdi"], 2),
  Smoke = rep(mean_values["Smoke"], 2),
  Drink = rep(mean_values["Drink"], 2),
  Mets = rep(mean_values["Mets"], 2),
  Dietscore = rep(mean_values["Dietscore"], 2),
  Sleeptime = rep(mean_values["Sleeptime"], 2),
  BMI = rep(mean_values["BMI"], 2),
  T2dhis = rep(mean_values["T2dhis"], 2),
  Cvdhis = rep(mean_values["Cvdhis"], 2),
  hiscancer = rep(mean_values["hiscancer"], 2),
  cancer_treat = rep(mean_values["cancer_treat"], 2),
  CMD_treat = rep(mean_values["CMD_treat"], 2),
  cancer_screen = rep(mean_values["cancer_screen"], 2),
  trans = 1:2
)
# Assign the transition matrix to the new dataset
attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
covs_P1 <- c("Age", "Sex", "Ethnicity", "Edu", 
             "Employed", "Tdi", "Smoke", "Drink", 
             "Mets", "Dietscore", "Sleeptime", "BMI",
             "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
             "CMD_treat", "cancer_screen")
# Expand the covariates for the new dataset
newd <- expand.covs(newd, covs_P1, longnames = FALSE)
# Define the strata for the transitions
newd$strata<-1:2
# Compute the multi-state probabilities
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
# Plot the survival curves
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)
# Estimate the mean values of probability and standard errors at 5, 10, 15 years
prob8_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob8_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob8_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
se8_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se8_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se8_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)

### 5. CA-CMD multimorbidity prognosis-------------------------------
track<-read.csv("CMD_CA.csv")
covariate<-read.csv("covariates_impute.csv")
data<-merge(track,covariate,by="n_eid")
# CA-CMD individuals
data <-data %>% filter (trans6==1)
# Define non-death outcome
data$trans11[data$trans9==1]=0
data$trans11[data$trans9==0]=1
# Select relevant columns from the dataset for each transition pattern
datatrans9<-data[,c(1,26,27,28,30:47)]
datatrans11<-data[,c(1,58,27,28,30:47)]
# Assign transition patterns to each dataset
datatrans9$transpattern<-1
datatrans11$transpattern<-2
# Rename the columns for each dataset
names(datatrans9)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
names(datatrans11)<-c("n_eid","status","timef","timeto","Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", "cancer_treat", "CMD_treat", "cancer_screen","failcode")
dataall1<-rbind(datatrans9,datatrans11)
# Calculate the time to the next event
dataall1$timeto<-dataall1$timeto-dataall1$timef
dataall1$timef<-0
library("mstate")
# Define the transition matrix for multi-state analysis (CA-CMD -> Death, non-Death)
tmat <- transMat(x = list(c(2,3), c(), c()),
                 names = c( "CA-CMD","Death","Non-Death"))
# Define the covariates for the model
covs_P <- c("Age", "Sex", "Ethnicity", "Edu", 
            "Employed", "Tdi", "Smoke", "Drink", 
            "Mets", "Dietscore", "Sleeptime", "BMI",
            "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
            "CMD_treat", "cancer_screen")
# Expand the covariates for the dataset
msebmt_allname <- expand.covs(dataall1,  covs_P , append = TRUE,longnames = FALSE)
# Get the column names of the expanded dataset
variables <- colnames(msebmt_allname[,24:59])
# Construct the formula 
formula_str <- paste("Surv(timef, timeto, status) ~", paste(variables, collapse = " + "), "+ strata(failcode)")
formula <- as.formula(formula_str)
the.expression_allname <- coxph(formula, data = msebmt_allname, method = 'breslow')
# Create a new dataset using the mean values of covariates
newd <- data.frame(
  Age = rep(mean_values["Age"], 2),
  Sex = rep(mean_values["Sex"], 2),
  Ethnicity = rep(mean_values["Ethnicity"], 2),
  Edu = rep(mean_values["Edu"], 2),
  Employed = rep(mean_values["Employed"], 2),
  Tdi = rep(mean_values["Tdi"], 2),
  Smoke = rep(mean_values["Smoke"], 2),
  Drink = rep(mean_values["Drink"], 2),
  Mets = rep(mean_values["Mets"], 2),
  Dietscore = rep(mean_values["Dietscore"], 2),
  Sleeptime = rep(mean_values["Sleeptime"], 2),
  BMI = rep(mean_values["BMI"], 2),
  T2dhis = rep(mean_values["T2dhis"], 2),
  Cvdhis = rep(mean_values["Cvdhis"], 2),
  hiscancer = rep(mean_values["hiscancer"], 2),
  cancer_treat = rep(mean_values["cancer_treat"], 2),
  CMD_treat = rep(mean_values["CMD_treat"], 2),
  cancer_screen = rep(mean_values["cancer_screen"], 2),
  trans = 1:2
)
# Assign the transition matrix to the new dataset
attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
covs_P1 <- c("Age", "Sex", "Ethnicity", "Edu", 
             "Employed", "Tdi", "Smoke", "Drink", 
             "Mets", "Dietscore", "Sleeptime", "BMI",
             "T2dhis", "Cvdhis", "hiscancer", "cancer_treat",
             "CMD_treat", "cancer_screen")
# Expand the covariates for the new dataset
newd <- expand.covs(newd, covs_P1, longnames = FALSE)
# Define the strata for the transitions
newd$strata<-1:2
# Compute the multi-state probabilities
msf.WW<-msfit(object = the.expression_allname ,newdata = newd,trans=tmat)
# Plot the survival curves
plot(msf.WW)
pt <- probtrans(msf.WW, predt = 0,)
pc<-summary(pt,from=1)
# Estimate the mean values of probability and standard errors at 5, 10, 15 years
prob9_5<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$upper2),3),")")
prob9_10<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$upper2),3),")")
prob9_15<-paste0(round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$pstate2),3)," (",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$lower2),3),", ",round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$upper2),3),")")
se9_5<-round(mean(pc[[1]][round(pc[[1]]$time,0)==5,]$se2),3)
se9_10<-round(mean(pc[[1]][round(pc[[1]]$time,0)==10,]$se2),3)
se9_15<-round(mean(pc[[1]][round(pc[[1]]$time,0)==15,]$se2),3)

proball<-rbind(prob1_5,prob1_10,prob1_15,prob2_5,prob2_10,prob2_15,prob3_5,prob3_10,prob3_15,
               prob4_5,prob4_10,prob4_15,prob5_5,prob5_10,prob5_15,prob6_5,prob6_10,prob6_15,
               prob7_5,prob7_10,prob7_15,prob8_5,prob8_10,prob8_15,prob9_5,prob9_10,prob9_15,
               se1_5,se1_10,se1_15,se2_5,se2_10,se2_15,
               se3_5,se3_10,se3_15,se4_5,se4_10,se4_15,
               se5_5,se5_10,se5_15,se6_5,se6_10,se6_15,
               se7_5,se7_10,se7_15,se8_5,se8_10,se8_15,
               se9_5,se9_10,se9_15
)

write.csv(proball,"all_trans.csv")