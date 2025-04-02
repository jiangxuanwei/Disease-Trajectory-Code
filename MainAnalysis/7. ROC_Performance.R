# ============================================================
# ROC Curves
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

library(riskRegression)  # For risk prediction evaluation
library(survival)        # For survival analysis
library(prodlim)         # For product-limit estimation
library(pec)             # For prediction error curves

# Load and prepare datasets ---------------------------------------------------

# Load pre-calculated risk scores for different biomarker types
data1 <- read.csv("PRSScore.csv")          # Polygenic Risk Scores (PRS)
data2 <- read.csv("MetaScore.csv")         # Metabolomics scores
data2 <- data2[,c(2,39:47)]               # Select relevant columns (ID + 9 transition scores)
data3 <- read.csv("ProtScore.csv")         # Proteomics scores  
data3 <- data3[,c(2,39:47)]               # Select relevant columns
corv <- read.csv("covariates_impute.csv")  # Covariates data

# Merge all datasets by participant ID
data <- merge(data1, data2, by="n_eid")
data <- merge(data, data3, by="n_eid")
datamerge <- merge(data, corv, by="n_eid")

# Calculate non-HDL cholesterol
datamerge$nonhdl <- datamerge$cho - datamerge$hdl

# Calculate median event times for each transition ----------------------------
# (Helps understand the time frame of each outcome)
median(datamerge[datamerge$outcome1==1,]$time1to)  # Transition 1: Healthy → CMD
median(datamerge[datamerge$outcome2==1,]$time2to)  # Transition 2: Healthy → Cancer
median(datamerge[datamerge$outcome3==1,]$time3to)  # Transition 3: Healthy → Death
median(datamerge[datamerge$outcome4==1,]$time4to)  # Transition 4: CMD → CMD+Cancer
median(datamerge[datamerge$outcome5==1,]$time5to)  # Transition 5: CMD → Death
median(datamerge[datamerge$outcome6==1,]$time6to)  # Transition 6: Cancer → CMD+Cancer
median(datamerge[datamerge$outcome7==1,]$time7to)  # Transition 7: Cancer → Death
median(datamerge[datamerge$outcome8==1,]$time8to)  # Transition 8: CMD+Cancer → Death
median(datamerge[datamerge$outcome9==1,]$time9to)  # Transition 9: CMD → Death

# Initialize result storage
resultAUC <- c()  # For AUC values
resultp <- c()    # For p-values of model comparisons

# ROC Analysis for each transition -------------------------------------------
for(i in 1:9){
  # Create PDF output for ROC curves (one per transition)
  pdf(paste0('outcomes 10',i,'.pdf'), 
      onefile = FALSE)
  
  # Define models for comparison:
  
  # 1. Base model (minimal covariates)
  Base <- paste0('fs1 <- coxph(Surv(time',i,'to,outcome',i,')~ Age+Sex, data=datamerge, x=TRUE)')
  
  # 2. Lifestyle model (base + lifestyle factors)
  Base1 <- paste0('fs2 <- coxph(Surv(time',i,'to,outcome',i,')~ Age+Sex+Smoke+Sleeptime+Mets+Dietscore+cancer_screen, data=datamerge, x=TRUE)')
  
  # 3. PRS model (genetic risk score only)
  PRS <- paste0('fs3 <- coxph(Surv(time',i,'to,outcome',i,')~ PRSscore',i,', data=datamerge, x=TRUE)')
  
  # 4. Clinical model (traditional clinical biomarkers)
  Chem <- paste0('fs4 <- coxph(Surv(time',i,'to,outcome',i,')~ Age+Sex+Smoke+sbp+nonhdl+glucose, data=datamerge, x=TRUE)')
  
  # 5. Metabolomics model
  Meta <- paste0('fs5 <- coxph(Surv(time',i,'to,outcome',i,')~ Metascore',i,', data=datamerge, x=TRUE)')  
  
  # 6. Proteomics model
  Prot <- paste0('fs6 <- coxph(Surv(time',i,'to,outcome',i,')~ Protscore',i,', data=datamerge, x=TRUE)')
  
  # 7. Combined model (all scores + covariates)
  Combined <- paste0('fs7 <- coxph(Surv(time',i,'to,outcome',i,')~ Metascore',i,'+Protscore',i,'+PRSscore',i,'+Smoke+Sleeptime+Mets+Dietscore+cancer_screen+sbp+nonhdl+glucose, data=datamerge, x=TRUE)')
  
  # Fit all models
  eval(parse(text = Base)) 
  eval(parse(text = Base1)) 
  eval(parse(text = Meta))
  eval(parse(text = Prot))
  eval(parse(text = PRS)) 
  eval(parse(text = Chem))
  eval(parse(text = Combined))
  
  # Evaluate models using time-dependent ROC at 10 years
  the.expression_allname <- paste0("xs<- Score(list(Base_model=fs1,Lifestyle_model=fs2,Clinical_model=fs4,PRS_Score=fs3,Met_Score=fs5,Prot_Score=fs6,Combined_Score=fs7), Hist(time",i,"to,outcome",i,")~1, data=datamerge, times=10, plots='roc', metrics='auc')")
  eval(parse(text = the.expression_allname))
  
  # Store results
  resultAUC <- cbind(resultAUC, xs[["AUC"]][["score"]][["AUC"]])  # AUC values
  resultp <- cbind(resultp, xs[["AUC"]][["contrasts"]][["p"]][1:11])  # p-values for comparisons
  
  # Generate ROC plot with custom colors
  p <- plotROC(xs, col = c('#7F7F80','#55719D',"#D6B0FF",'#5C996C',"#E2C7BA",'#FF9045','#BF1E20'))
  dev.off()  # Close PDF device
}

# Save results --------------------------------------------------------------
write.csv(resultAUC, "AUC 10y.csv")
write.csv(resultp, "AUC-p 10y.csv")

### AUC in 15 years follow the same analysis pipeline as described above