README

This repository contains code used to generate the results of the article titled "Disease Trajectories of Cardiometabolic Diseases and Cancer: Transition Patterns, Multiomics Signatures, Prognosis and Prediction".

Codes for producing publication figures are provided for transparency and in the spirit of scientific collaboration. We will not be able to answer questions about the details of these codes.

Below is an overview of the contents of each folder.

1. DataPreparation Folder

   This folder includes all the main data preprocessing steps and core data file generation scripts. The scripts in this folder are essential for cleaning, imputing, and defining variables, as well as creating the necessary data files for further analysis.

   Files in DataPreparation:
   
   Disease_ICD10_Code.R
      - This script defines the primary outcome variables and records the occurrence time of the diseases using ICD-10 codes.
   
   Impute_Corv_Data.R
      - This script handles the imputation process for missing covariate data.
   
   UKB_Disease_Trajectory_Analysis.R
      - This script defines the transition outcomes and the start and end times for each transition event in the analysis.
   
   Subgroup_Definition.R
      - This script defines and generates various subgroups based on specific criteria, preparing the dataset for subgroup analysis.

2. MainAnalysis Folder

   This folder includes all the scripts used for performing the main analyses and generating the primary results. The analysis includes the calculation of transition probabilities, survival analysis, multiomics signature identification, and other essential steps in the research.

   Files in MainAnalysis:
   
   Transition_Probability_Calculation.R
      - This script calculates the transition probabilities and standard errors for three distinct stages of disease progression.
   
   KM_Plot and RMST_Calculation.R
      - This script generates Kaplan-Meier survival curves and calculates Restricted Mean Survival Time (RMST) for different stages.
   
   Multiomics_Signatures_Identification.R
      - This script identifies important multiomics signatures for each transition stage using a multistate model.
   
   KEGG_GO.R
      - This script performs KEGG and Gene Ontology (GO) enrichment analysis on the identified multiomics signatures to understand their biological significance.
   
   Principal_Component.R
      - This script performs Principal Component Analysis (PCA) based on the major multiomics signatures to reduce dimensionality and visualize the data.
   
   LASSO_Model_For_Prediction.R
      - This script implements the LASSO regression model to identify predictive features (predictors) and evaluates the prediction performance of the model.
   
   ROC_Performance.R
      - This script generates Receiver Operating Characteristic (ROC) curves and calculates the Area Under the Curve (AUC) to assess the prediction model's performance.

3. Figures and Tables Folder

   This folder contains the scripts used to generate the key figures and tables for the analysis. The visualizations and tabular summaries are critical for presenting the findings in the paper.

License
This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa]. [![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]
