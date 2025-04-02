# ============================================================
# Sankey Diagram
# Author: JXW
# Last Modified: Feb 2025
# ============================================================

library(networkD3)
links <- read.csv("links_path1.csv",header = TRUE,stringsAsFactors = FALSE)
nodes <- read.csv("nodes_path1.csv",header = TRUE,stringsAsFactors = FALSE)
# Prepare source and target IDs for the Sankey diagram
# Match source names to node names and convert to zero-based index (required by D3)
links$IDsource <- match(links$source, nodes$nodes)-1 
# Same for target names
links$IDtarget <- match(links$target, nodes$nodes)-1
# Define color scale using D3.js syntax
my_color <- 'd3.scaleOrdinal() .domain(["Health","CMD-1","CA-1","CMD-2","CA-2"]) .range(["#3864A4","#9F1F24","#9B67A9","#9F1F24","#9B67A9"])'
# Create node groups for coloring purposes
nodes$group[nodes$Name %in% c("Health")]<-"A"
nodes$group[nodes$Name %in% c("CMD-1","CMD-2","Only CMD","Only CA-CMD")]<-"B"
nodes$group[nodes$Name %in% c("CA-1","CA-2","Only CA","Only CMD-CA")]<-"C"
library(networkD3)
sankeyNetwork(Links = links, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "Value", NodeID = "nodes", fontSize = 8,
                       colourScale="my_color",
                       NodeGroup="group",
                       LinkGroup = 'source', 
                       sinksRight=FALSE)