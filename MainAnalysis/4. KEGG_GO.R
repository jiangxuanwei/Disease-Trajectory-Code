# ============================================================
# KEGG and GO pathways
# Author: YGR
# Last Modified: Feb 2025
# ============================================================

# Load necessary libraries for enrichment analysis and Excel file handling
library(clusterProfiler)  # For enrichment analysis
library(org.Hs.eg.db)     # For mapping protein IDs to Entrez IDs
library(openxlsx)         # For creating and writing Excel files

# Read in the significant protein data from CSV file
protein <- read.csv("Prot_significant.csv")

# Create a new workbook to store the results
wb <- createWorkbook()

# Add a summary sheet to the workbook
addWorksheet(wb, sheetName = "Summary")

# Initialize an empty dataframe to store the summary of enrichment results
summary_data <- data.frame(
  Variable = character(),
  Input_Proteins = integer(),
  Mapped_Proteins = integer(),
  GO_Pathways = integer(),
  KEGG_Pathways = integer(),
  stringsAsFactors = FALSE
)

# Define the variables to analyze (protein lists in the dataset)
variables <- c("id1", "id2", "id3", "id4", "id5", "id6", "id7", "id8", "id9")

# Loop over each variable to perform enrichment analysis
for (variable in variables) {
  
  # Check if the variable exists in the protein dataset
  if (!(variable %in% names(protein))) {
    cat("Warning: Variable", variable, "does not exist in the dataset, skipping.\n")
    next  # Skip this iteration if the variable is missing
  }
  
  # Extract the protein list for the current variable
  protein_list <- protein[[variable]]
  
  # If the protein list has no valid entries, skip this iteration
  if (length(protein_list[!is.na(protein_list) & protein_list != ""]) == 0) {
    cat("Warning: Variable", variable, "has no valid protein data, skipping.\n")
    next  # Skip if no valid proteins are found
  }
  
  # Perform the enrichment analysis using the protein list
  results <- enrich_analysis(protein_list, variable)
  
  # Count the number of input proteins (non-missing and non-empty entries)
  input_count <- length(protein_list[!is.na(protein_list) & protein_list != ""])
  
  # Map the protein list to Entrez IDs and count the number of mapped proteins
  mapped_count <- if (!is.null(results$GO) || !is.null(results$KEGG)) {
    nrow(bitr(toupper(protein_list[!is.na(protein_list) & protein_list != ""]), 
              fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  } else {
    0  # If no enrichment results, set mapped_count to 0
  }
  
  # Count the number of pathways identified for Gene Ontology (GO) and KEGG
  go_count <- if (!is.null(results$GO)) nrow(results$GO) else 0
  kegg_count <- if (!is.null(results$KEGG)) nrow(results$KEGG) else 0
  
  # Add a row to the summary data for the current variable
  summary_data <- rbind(summary_data, data.frame(
    Variable = variable,
    Input_Proteins = input_count,
    Mapped_Proteins = mapped_count,
    GO_Pathways = go_count,
    KEGG_Pathways = kegg_count,
    stringsAsFactors = FALSE
  ))
  
  # If GO pathways are found, add a new worksheet with the GO results
  if (!is.null(results$GO) && nrow(results$GO) > 0) {
    sheet_name <- paste0(variable, "_GO")  # Create sheet name based on variable
    addWorksheet(wb, sheetName = sheet_name)  # Add the sheet to the workbook
    writeData(wb, sheet = sheet_name, x = results$GO)  # Write GO results to the sheet
    setColWidths(wb, sheet = sheet_name, cols = 1:ncol(results$GO), widths = "auto")  # Auto-adjust column widths
    setColWidths(wb, sheet = sheet_name, cols = which(colnames(results$GO) == "Proteins"), widths = 50)  # Set wider column for "Proteins"
  }
  
  # If KEGG pathways are found, add a new worksheet with the KEGG results
  if (!is.null(results$KEGG) && nrow(results$KEGG) > 0) {
    sheet_name <- paste0(variable, "_KEGG")  # Create sheet name based on variable
    addWorksheet(wb, sheetName = sheet_name)  # Add the sheet to the workbook
    writeData(wb, sheet = sheet_name, x = results$KEGG)  # Write KEGG results to the sheet
    setColWidths(wb, sheet = sheet_name, cols = 1:ncol(results$KEGG), widths = "auto")  # Auto-adjust column widths
    setColWidths(wb, sheet = sheet_name, cols = which(colnames(results$KEGG) == "Proteins"), widths = 50)  # Set wider column for "Proteins"
  }
}

# Write the summary data to the "Summary" sheet in the workbook
writeData(wb, sheet = "Summary", x = summary_data)

# Auto-adjust column widths in the summary sheet
setColWidths(wb, sheet = "Summary", cols = 1:ncol(summary_data), widths = "auto")

# Save the workbook to a file
saveWorkbook(wb, "Protein_Enrichment_Results1.xlsx", overwrite = TRUE)
