# ============================================================
# Venn Plot
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

# Install and load VennDiagram package
library(VennDiagram)

# Read data
df <- read.csv("D:/科研数据/课题数据/UKB分析/疾病轨迹--本地/1. 修回文件/8.Markers/birthplace/Prot_significant.csv")

# Extract unique IDs from each column, removing NA values
set1 <- unique(na.omit(df$id1))
set2 <- unique(na.omit(df$id2))
set3 <- unique(na.omit(df$id3))
set4 <- unique(na.omit(df$id4))
set5 <- unique(na.omit(df$id5))
set6 <- unique(na.omit(df$id6))
set7 <- unique(na.omit(df$id7))
set8 <- unique(na.omit(df$id8))
set9 <- unique(na.omit(df$id9))

# Create first Venn diagram (3 sets) and save to PDF
pdf("D:/科研数据/课题数据/UKB分析/疾病轨迹--本地/1. 修回文件/8.Markers/venn_diagram_baseline.pdf", width = 10, height = 10)
sets_list <- list(
  "Set1" = set1,
  "Set2" = set2,
  "Set3" = set3
)

# Generate Venn diagram
venn.plot <- venn.diagram(
  x = sets_list,
  category.names = names(sets_list),
  col = c("#BF1E20", "#FF9045", "#7F7F80"),  # Border colors
  fill = c("#DD858C", "#F2BD91", "#CECECE"), # Fill colors
  filename = NULL,  # Don't save to file directly
  output = TRUE     # Return grob object
)

# Draw the plot
grid.draw(venn.plot)
dev.off()  # Close PDF device

# Create second Venn diagram (4 sets) and save to PDF
pdf("D:/科研数据/课题数据/UKB分析/疾病轨迹--本地/1. 修回文件/8.Markers/birthplace/venn_diagram_1.pdf", width = 10, height = 10)
sets_list <- list(
  "Set4" = set4,
  "Set5" = set5,
  "Set6" = set6,
  "Set7" = set7
)

# Generate Venn diagram
venn.plot <- venn.diagram(
  x = sets_list,
  category.names = names(sets_list),
  col = c("#BF1E20", "#FF9045", "#DB7BFF", "#6155FF"),
  fill = c("#DD858C", "#F2BD91", "#EFCAFF", "#CECAFF"),
  filename = NULL,
  output = TRUE
)

# Draw the plot
grid.draw(venn.plot)
dev.off()  # Close PDF device

# Create third Venn diagram (2 sets) and save to PDF
pdf("D:/科研数据/课题数据/UKB分析/疾病轨迹--本地/1. 修回文件/8.Markers/birthplace/venn_diagram_2.pdf", width = 10, height = 10)
sets_list <- list(
  "Set8" = set8,
  "Set9" = set9
)

# Generate Venn diagram
venn.plot <- venn.diagram(
  x = sets_list,
  category.names = names(sets_list),
  col = c("#DB7BFF", "#6155FF"),
  fill = c("#EFCAFF", "#CECAFF"),
  filename = NULL,
  output = TRUE
)

# Draw the plot
grid.draw(venn.plot)
dev.off()  # Close PDF device