# ============================================================
# ICD10-Code
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

### 1. Non-cancer definition-------------------------------------------
# T2D
data$time_dia=NA
for(i in 1:nrow(data)){
  TDtime<-data.frame()
  for(j in 0:258){
    if(data[, paste0("s_41270_0_", j)][i] %in% c("E110","E111","E112","E113","E114","E115","E116","E117","E118","E119")) {
      data$event_dia[i] <- 1
      TDtime<-rbind (TDtime,data[, paste0("s_41280_0_", j)][i])
      data$time_dia[i] <- (min(TDtime) - data$s_53_0_0[i])/365
    }
  }
}

# CVD
CVD_codes <- c("I20", "I21", "I22", "I23", "I24", "I25", "I50", "I60", "I61", "I62", "I63", "I64")
data$event_CVD=0
data$time_CVD=NA
for(i in 1:nrow(data)){
  CVDtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", CVD_codes, "$|^", CVD_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_CVD[i] <- 1
      CVDtime<-rbind(CVDtime,data[, paste0("s_41280_0_", j)][i])
      data$time_CVD[i] <- (min(CVDtime) - data$s_53_0_0[i])/365
    }
  }
}

# CAD
CAD_codes <- c("I20", "I21", "I22", "I23", "I24", "I25")
data$event_CAD=0
data$time_CAD=NA
for(i in 1:nrow(data)){
  CADtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", CAD_codes, "$|^", CAD_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_CAD[i] <- 1
      CADtime<-rbind(CADtime,data[, paste0("s_41280_0_", j)][i])
      data$time_CAD[i] <- (min(CADtime) - data$s_53_0_0[i])/365
    }
  }
}

# HF
HF_codes <- c("I50")
data$event_HF=0
data$time_HF=NA
for(i in 1:nrow(data)){
  HFtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", HF_codes, "$|^", HF_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_HF[i] <- 1
      HFtime<-rbind(HFtime,data[, paste0("s_41280_0_", j)][i])
      data$time_HF[i] <- (min(HFtime) - data$s_53_0_0[i])/365
    }
  }
}

# stroke
stroke_codes <- c("I60", "I61", "I62", "I63", "I64")
data$event_stroke=0
data$time_stroke=NA
for(i in 1:nrow(data)){
  stroketime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", stroke_codes, "$|^", stroke_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_stroke[i] <- 1
      stroketime<-rbind(stroketime,data[, paste0("s_41280_0_", j)][i])
      data$time_stroke[i] <- (min(stroketime) - data$s_53_0_0[i])/365
    }
  }
}

### 2. Cancer definition-------------------------------------------
# All Cancer
data$event_CA=0
data$time_CA=NA
for(i in 1:nrow(data)){
    if(!is.na(data$s_40005_0_0[i])) {
      data$event_CA[i] <- 1
      data$time_CA[i] <- (data$s_40005_0_0[i] - data$s_53_0_0[i])/365
  }
}

# lung Cancer
LungCA_codes <- c("C33", "C34")
data$event_LungCA=0
data$time_LungCA=NA
for(i in 1:nrow(data)){
  LungCAtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", LungCA_codes, "$|^", LungCA_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_LungCA[i] <- 1
      LungCAtime<-rbind(LungCAtime,data[, paste0("s_41280_0_", j)][i])
      data$time_LungCA[i] <- (min(LungCAtime) - data$s_53_0_0[i])/365
    }
  }
}

# Skin Cancer
SkinCA_codes <- c("C44")
data$event_SkinCA=0
data$time_SkinCA=NA
for(i in 1:nrow(data)){
  SkinCAtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", SkinCA_codes, "$|^", SkinCA_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_SkinCA[i] <- 1
      SkinCAtime<-rbind(SkinCAtime,data[, paste0("s_41280_0_", j)][i])
      data$time_SkinCA[i] <- (min(SkinCAtime) - data$s_53_0_0[i])/365
    }
  }
}

# Prostate Cancer
ProsCA_codes <- c("C61")
data$event_ProsCA=0
data$time_ProsCA=NA
for(i in 1:nrow(data)){
  ProsCAtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", ProsCA_codes, "$|^", ProsCA_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_ProsCA[i] <- 1
      ProsCAtime<-rbind(ProsCAtime,data[, paste0("s_41280_0_", j)][i])
      data$time_ProsCA[i] <- (min(ProsCAtime) - data$s_53_0_0[i])/365
    }
  }
}

# Breast Cancer
BreastCA_codes <- c("C50")
data$event_BreastCA=0
data$time_BreastCA=NA
for(i in 1:nrow(data)){
  BreastCAtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", BreastCA_codes, "$|^", BreastCA_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_BreastCA[i] <- 1
      BreastCAtime<-rbind(BreastCAtime,data[, paste0("s_41280_0_", j)][i])
      data$time_BreastCA[i] <- (min(BreastCAtime) - data$s_53_0_0[i])/365
    }
  }
}

# Stomtate Cancer
StomCA_codes <- c("C16")
data$event_StomCA=0
data$time_StomCA=NA
for(i in 1:nrow(data)){
  StomCAtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", StomCA_codes, "$|^", StomCA_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_StomCA[i] <- 1
      StomCAtime<-rbind(StomCAtime,data[, paste0("s_41280_0_", j)][i])
      data$time_StomCA[i] <- (min(StomCAtime) - data$s_53_0_0[i])/365
    }
  }
}

# Liver Cancer
data$event_LiverCA=0
data$time_LiverCA=NA
for(i in 1:nrow(data)){
  LiverCAtime<-data.frame()
  for(j in 0:258){
    if(data[, paste0("s_41270_0_", j)][i] %in% c('K760','C220',"C221","C222","C223","C224","C227","C229")) {
      data$event_LiverCA[i] <- 1
      LiverCAtime<-rbind(LiverCAtime,data[, paste0("s_41280_0_", j)][i])
      data$time_LiverCA[i] <- (min(LiverCAtime) - data$s_53_0_0[i])/365
    }
  }
}

# Esophageal Cancer
EsophagealCA_codes <- c("C15")
data$event_EsophagealCA=0
data$time_EsophagealCA=NA
for(i in 1:nrow(data)){
  EsophagealCAtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", EsophagealCA_codes, "$|^", EsophagealCA_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_EsophagealCA[i] <- 1
      EsophagealCAtime<-rbind(EsophagealCAtime,data[, paste0("s_41280_0_", j)][i])
      data$time_EsophagealCA[i] <- (min(EsophagealCAtime) - data$s_53_0_0[i])/365
    }
  }
}

# Pancreatic Cancer
PancreaticCA_codes <- c("C25")
data$event_PancreaticCA=0
data$time_PancreaticCA=NA
for(i in 1:nrow(data)){
  PancreaticCAtime<-data.frame()
  for(j in 0:258){
    if(any(grepl(paste0("(^", PancreaticCA_codes, "$|^", PancreaticCA_codes, "[0-9]$)", collapse = "|"),data[, paste0("s_41270_0_", j)][i]))) {
      data$event_PancreaticCA[i] <- 1
      PancreaticCAtime<-rbind(PancreaticCAtime,data[, paste0("s_41280_0_", j)][i])
      data$time_PancreaticCA[i] <- (min(PancreaticCAtime) - data$s_53_0_0[i])/365
    }
  }
}

### 3. All-cause mortality-------------------------------------------
for(i in 1:nrow(data)){
    if (!is.na(data$s_40000_0_0[i]))
    {data$ death2_icd10[i]=1
    data$ death2_time[i]=(data$s_40000_0_0[i]-data$s_53_0_0[i])/365
    }
    else if (is.na(data$s_40000_0_0[i]))
    {data$ death2_icd10[i]=0
    data$ death2_time[i]=(23303-data$s_53_0_0[i])/365
    }
}

