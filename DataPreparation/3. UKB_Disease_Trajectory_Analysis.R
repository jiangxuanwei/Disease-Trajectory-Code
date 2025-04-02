# ============================================================
# Disease Trajectory Analysis Script
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

# Load dplyr package for data manipulation
library(dplyr)

# Read data from CSV file
data<-read.csv("D:/track.csv")

# Create composite cardiovascular disease (CD) variable
# CD = 1 if either diabetes (event_dia) or CVD (event_CVD) occurred
data$event_CD=0
data$event_CD[data$event_dia==1|data$event_CVD==1]=1

# Calculate time to first CD event (minimum of diabetes or CVD time)
data$time_CD=pmin(data$time_dia,data$time_CVD,na.rm=T)

# Filter data:
# 1. Keep records where time_CD is NA or >0
# 2. Keep records where time_CA is NA or >0
# 3. Exclude records where time_CD equals time_CA (unless NA)
data1<-data %>%filter(is.na(time_CD)|time_CD>0) %>%filter(is.na(time_CA)|time_CA>0) %>%
  filter(time_CD-time_CA!=0|is.na(time_CD-time_CA)) 

# Add small constant (0.1) to death time to avoid numerical issues
data1$death2_time<- data1$death2_time+0.1

# Define transition states for multi-state model
# Loops through each subject to determine their disease trajectory
for(i in 1:nrow(data1)){
  # Case 1: No CD, No CA, No death (censored)
  if (data1$event_CD[i]==0 & data1$event_CA[i]==0 & data1$death2_icd10[i]==0 )
  {data1$ trans1[i]=0
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$death2_time[i]
  data1$ trans2[i]=0
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$death2_time[i]
  data1$ trans3[i]=0
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$death2_time[i]  
  data1$ trans4[i]=0
  data1$ time4from[i]=0
  data1$ time4to[i]=data1$death2_time[i]  
  data1$ trans5[i]=0
  data1$ time5from[i]=0
  data1$ time5to[i]=data1$death2_time[i] 
  data1$ trans6[i]=0
  data1$ time6from[i]=0
  data1$ time6to[i]=data1$death2_time[i]  
  data1$ trans7[i]=0
  data1$ time7from[i]=0
  data1$ time7to[i]=data1$death2_time[i]
  data1$ trans8[i]=0
  data1$ time8from[i]=0
  data1$ time8to[i]=data1$death2_time[i]  
  data1$ trans9[i]=0
  data1$ time9from[i]=0
  data1$ time9to[i]=data1$death2_time[i]
  }
  
  # Case 2: No CD, No CA, but death occurred
  if (data1$event_CD[i]==0 & data1$event_CA[i]==0 & data1$death2_icd10[i]==1)
  {data1$ trans1[i]=0
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$death2_time[i]
  data1$ trans2[i]=0
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$death2_time[i]
  data1$ trans3[i]=1
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$death2_time[i]  
  data1$ trans4[i]=0
  data1$ time4from[i]=0
  data1$ time4to[i]=data1$death2_time[i]  
  data1$ trans5[i]=0
  data1$ time5from[i]=0
  data1$ time5to[i]=data1$death2_time[i] 
  data1$ trans6[i]=0
  data1$ time6from[i]=0
  data1$ time6to[i]=data1$death2_time[i]  
  data1$ trans7[i]=0
  data1$ time7from[i]=0
  data1$ time7to[i]=data1$death2_time[i]
  data1$ trans8[i]=0
  data1$ time8from[i]=0
  data1$ time8to[i]=data1$death2_time[i]  
  data1$ trans9[i]=0
  data1$ time9from[i]=0
  data1$ time9to[i]=data1$death2_time[i]
  }
  
  # Case 3: No CD, CA occurred, No death
  if (data1$event_CD[i]==0 & data1$event_CA[i]==1 & data1$death2_icd10[i]==0)
  {data1$ trans1[i]=0
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$time_CA[i]
  data1$ trans2[i]=1
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$time_CA[i]
  data1$ trans3[i]=0
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$time_CA[i]
  data1$ trans4[i]=0
  data1$ time4from[i]=0
  data1$ time4to[i]=data1$death2_time[i]  
  data1$ trans5[i]=0
  data1$ time5from[i]=0
  data1$ time5to[i]=data1$death2_time[i] 
  data1$ trans6[i]=0
  data1$ time6from[i]=data1$time_CA[i]
  data1$ time6to[i]=data1$death2_time[i]  
  data1$ trans7[i]=0
  data1$ time7from[i]=data1$time_CA[i]
  data1$ time7to[i]=data1$death2_time[i] 
  data1$ trans8[i]=0
  data1$ time8from[i]=0
  data1$ time8to[i]=data1$death2_time[i]
  data1$ trans9[i]=0
  data1$ time9from[i]=data1$time_CA[i]
  data1$ time9to[i]=data1$death2_time[i] 
  }
  
  # Case 4: CD occurred, No CA, No death
  if (data1$event_CD[i]==1 & data1$event_CA[i]==0 & data1$death2_icd10[i]==0)
  {data1$ trans1[i]=1
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$time_CD[i]
  data1$ trans2[i]=0
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$time_CD[i]
  data1$ trans3[i]=0
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$time_CD[i]
  data1$ trans4[i]=0
  data1$ time4from[i]=data1$time_CD[i]
  data1$ time4to[i]=data1$death2_time[i]  
  data1$ trans5[i]=0
  data1$ time5from[i]=data1$time_CD[i]
  data1$ time5to[i]=data1$death2_time[i] 
  data1$ trans6[i]=0
  data1$ time6from[i]=0
  data1$ time6to[i]=data1$death2_time[i]  
  data1$ trans7[i]=0
  data1$ time7from[i]=0
  data1$ time7to[i]=data1$death2_time[i] 
  data1$ trans8[i]=0
  data1$ time8from[i]=data1$time_CD[i]
  data1$ time8to[i]=data1$death2_time[i]
  data1$ trans9[i]=0
  data1$ time9from[i]=0
  data1$ time9to[i]=data1$death2_time[i] 
  }
  
  # Case 5: CD occurred first, then CA, No death
  if (data1$event_CD[i]==1 & data1$event_CA[i]==1 & data1$death2_icd10[i]==0 & data1$time_CA[i]>data1$time_CD[i])
  {data1$ trans1[i]=1
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$time_CD[i]
  data1$ trans2[i]=0
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$time_CD[i]
  data1$ trans3[i]=0
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$time_CD[i]
  data1$ trans4[i]=1
  data1$ time4from[i]=data1$time_CD[i]
  data1$ time4to[i]=data1$time_CA[i]
  data1$ trans5[i]=0
  data1$ time5from[i]=data1$time_CD[i]
  data1$ time5to[i]=data1$time_CA[i] 
  data1$ trans6[i]=0
  data1$ time6from[i]=0
  data1$ time6to[i]=data1$death2_time[i]  
  data1$ trans7[i]=0
  data1$ time7from[i]=0
  data1$ time7to[i]=data1$death2_time[i] 
  data1$ trans8[i]=0
  data1$ time8from[i]=data1$time_CA[i]
  data1$ time8to[i]=data1$death2_time[i]
  data1$ trans9[i]=0
  data1$ time9from[i]=0
  data1$ time9to[i]=data1$death2_time[i] 
  }
  
  # Case 6: CA occurred first, then CD, No death
  if (data1$event_CD[i]==1 & data1$event_CA[i]==1 & data1$death2_icd10[i]==0 & data1$time_CA[i]<data1$time_CD[i])
  {data1$ trans1[i]=0
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$time_CA[i]
  data1$ trans2[i]=1
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$time_CA[i]
  data1$ trans3[i]=0
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$time_CA[i]
  data1$ trans4[i]=0
  data1$ time4from[i]=0
  data1$ time4to[i]=data1$death2_time[i] 
  data1$ trans5[i]=0
  data1$ time5from[i]=0
  data1$ time5to[i]=data1$death2_time[i] 
  data1$ trans6[i]=1
  data1$ time6from[i]=data1$time_CA[i]
  data1$ time6to[i]=data1$time_CD[i]  
  data1$ trans7[i]=0
  data1$ time7from[i]=data1$time_CA[i]
  data1$ time7to[i]=data1$time_CD[i] 
  data1$ trans8[i]=0
  data1$ time8from[i]=0
  data1$ time8to[i]=data1$death2_time[i]
  data1$ trans9[i]=0
  data1$ time9from[i]=data1$time_CD[i]
  data1$ time9to[i]=data1$death2_time[i] 
  }
  
  # Case 7: No CD, CA occurred, and death occurred
  if (data1$event_CD[i]==0 & data1$event_CA[i]==1 & data1$death2_icd10[i]==1 )
  {data1$ trans1[i]=0
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$time_CA[i]
  data1$ trans2[i]=1
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$time_CA[i]
  data1$ trans3[i]=0
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$time_CA[i]
  data1$ trans4[i]=0
  data1$ time4from[i]=0
  data1$ time4to[i]=data1$death2_time[i] 
  data1$ trans5[i]=0
  data1$ time5from[i]=0
  data1$ time5to[i]=data1$death2_time[i] 
  data1$ trans6[i]=0
  data1$ time6from[i]=data1$time_CA[i]
  data1$ time6to[i]=data1$death2_time[i]  
  data1$ trans7[i]=1
  data1$ time7from[i]=data1$time_CA[i]
  data1$ time7to[i]=data1$death2_time[i] 
  data1$ trans8[i]=0
  data1$ time8from[i]=0
  data1$ time8to[i]=data1$death2_time[i]
  data1$ trans9[i]=0
  data1$ time9from[i]=data1$time_CA[i]
  data1$ time9to[i]=data1$death2_time[i] 
  }
  
  # Case 8: CD occurred, No CA, but death occurred
  if (data1$event_CD[i]==1 & data1$event_CA[i]==0 & data1$death2_icd10[i]==1 )
  {data1$ trans1[i]=1
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$time_CD[i]
  data1$ trans2[i]=0
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$time_CD[i]
  data1$ trans3[i]=0
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$time_CD[i]
  data1$ trans4[i]=0
  data1$ time4from[i]=data1$time_CD[i]
  data1$ time4to[i]=data1$death2_time[i]
  data1$ trans5[i]=1
  data1$ time5from[i]=data1$time_CD[i]
  data1$ time5to[i]=data1$death2_time[i] 
  data1$ trans6[i]=0
  data1$ time6from[i]=0
  data1$ time6to[i]=data1$death2_time[i]  
  data1$ trans7[i]=0
  data1$ time7from[i]=0
  data1$ time7to[i]=data1$death2_time[i] 
  data1$ trans8[i]=0
  data1$ time8from[i]=data1$time_CD[i]
  data1$ time8to[i]=data1$death2_time[i]
  data1$ trans9[i]=0
  data1$ time9from[i]=0
  data1$ time9to[i]=data1$death2_time[i] 
  }
  
  # Case 9: CA occurred first, then CD, and death occurred
  if (data1$event_CD[i]==1 & data1$event_CA[i]==1 & data1$death2_icd10[i]==1 & data1$time_CA[i]<data1$time_CD[i])
  {data1$ trans1[i]=0
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$time_CA[i]
  data1$ trans2[i]=1
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$time_CA[i]
  data1$ trans3[i]=0
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$time_CA[i]
  data1$ trans4[i]=0
  data1$ time4from[i]=0
  data1$ time4to[i]=data1$death2_time[i] 
  data1$ trans5[i]=0
  data1$ time5from[i]=0
  data1$ time5to[i]=data1$death2_time[i] 
  data1$ trans6[i]=1
  data1$ time6from[i]=data1$time_CA[i]
  data1$ time6to[i]=data1$time_CD[i]  
  data1$ trans7[i]=0
  data1$ time7from[i]=data1$time_CA[i]
  data1$ time7to[i]=data1$time_CD[i] 
  data1$ trans8[i]=0
  data1$ time8from[i]=0
  data1$ time8to[i]=data1$death2_time[i]
  data1$ trans9[i]=1
  data1$ time9from[i]=data1$time_CD[i]
  data1$ time9to[i]=data1$death2_time[i] 
  }
  
  # Case 10: CD occurred first, then CA, and death occurred
  if (data1$event_CD[i]==1 & data1$event_CA[i]==1 & data1$death2_icd10[i]==1 & data1$time_CA[i]>data1$time_CD[i])
  {data1$ trans1[i]=1
  data1$ time1from[i]=0
  data1$ time1to[i]=data1$time_CD[i]
  data1$ trans2[i]=0
  data1$ time2from[i]=0
  data1$ time2to[i]=data1$time_CD[i]
  data1$ trans3[i]=0
  data1$ time3from[i]=0
  data1$ time3to[i]=data1$time_CD[i]
  data1$ trans4[i]=1
  data1$ time4from[i]=data1$time_CD[i]
  data1$ time4to[i]=data1$time_CA[i]
  data1$ trans5[i]=0
  data1$ time5from[i]=data1$time_CD[i]
  data1$ time5to[i]=data1$time_CA[i] 
  data1$ trans6[i]=0
  data1$ time6from[i]=0
  data1$ time6to[i]=data1$death2_time[i]  
  data1$ trans7[i]=0
  data1$ time7from[i]=0
  data1$ time7to[i]=data1$death2_time[i] 
  data1$ trans8[i]=1
  data1$ time8from[i]=data1$time_CA[i]
  data1$ time8to[i]=data1$death2_time[i]
  data1$ trans9[i]=0
  data1$ time9from[i]=0
  data1$ time9to[i]=data1$death2_time[i] 
  }
}

# Save processed data to CSV file
write.csv(data1,"D:/CMD_CA.csv")
