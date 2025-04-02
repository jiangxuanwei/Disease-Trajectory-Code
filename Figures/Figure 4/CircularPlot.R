# ============================================================
# Circular Plot
# Author: YGR
# Last Modified: Feb 2025
# ============================================================

library(openxlsx)
library(dplyr)
library(stringr)

# Load the workbook and read data
wb <- loadWorkbook("Pathway_go.xlsx")
data <- read.xlsx(wb, "Pathway")  # Read Pathway data sheet
abbre <- read.xlsx(wb, "Abbre")   # Read Abbreviation data sheet

# Merge the two datasets by 'Pathway'
data <- merge(data, abbre, by = "Pathway")

# Calculate -log10(p-value) and store in a new column 'value'
data <- data %>%
  mutate(value = -log10(as.numeric(Pvalue)))

# Define the groups for the data (transitions)
datagroup <- c("Transition1", "Transition2" ,"Transition3" ,"Transition4", "Transition5", "Transition6", "Transition7", "Transition8", "Transition9")

# Sort the data by group and value in descending order
data_sorted <- data %>%
  group_by(group) %>%
  arrange(desc(value)) %>%
  ungroup()

# Create empty placeholder data
empty_data <- tibble(
  'group' = datagroup,
  'individual' = paste0('empty_individual_', seq_along(datagroup)),
  'value' = 0
)

# Combine the original data with the empty placeholder data for each group
allplotdata <- NULL
for(g in datagroup) {
  # Get the empty placeholder data for the current group
  empty_g <- empty_data %>% filter(group == g)
  # Get the sorted data for the current group
  data_g <- data_sorted %>% filter(group == g)
  # Combine the empty placeholder and actual data for the group
  allplotdata <- bind_rows(allplotdata, empty_g, data_g)
}

# Add position information for plotting
allplotdata <- allplotdata %>% 
  mutate(xid = 1:n()) %>%  # Create unique IDs for each data point
  mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>%  # Calculate the angle for each point in the circular plot
  mutate(hjust = ifelse(angle < -90, 1, 0)) %>%  # Adjust horizontal alignment based on angle
  mutate(angle = ifelse(angle < -90, angle + 180, angle))  # Adjust angle if it's less than -90 degrees

# Extract and adjust empty placeholder data for segment labeling
firstxid <- which(str_detect(allplotdata$individual, pattern = "empty_individual")) 
segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to) / 2))

# Define custom y-axis positions for labels
coordy <- tibble('coordylocation' = seq(from = min(allplotdata$value), to = max(allplotdata$value), 10),
                 'coordytext' = as.character(round(coordylocation, 2)),
                 'x' = 1)

# Define grid data for custom axis gridlines
griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)

# Create the circular plot using ggplot2
p <- ggplot() + 
  geom_bar(data = allplotdata, aes(x = xid, y = value, fill = group), stat = 'identity') +  # Bar plot for each data point
  geom_text(data = allplotdata %>% filter(!str_detect(individual, pattern = "empty_individual")), 
            aes(x = xid, label = individual, y = value + 3, angle = angle, hjust = hjust),
            color = "black", fontface = "bold", alpha = 0.6, size = 2.5) +  # Labels for each data point
  geom_segment(data = segment_data, aes(x = from, xend = to), y = -5, yend = -5) +  # Segments for the transitions
  geom_text(data = segment_data, aes(x = labelx, label = label), y = -15) +  # Labels for each transition segment
  geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
            color = "grey", size = 3, angle = 0, fontface = "bold") +  # y-axis labels
  geom_segment(data = griddata, 
               aes(x = locationx - 0.5, xend = locationx + 0.5, y = locationy, yend = locationy),
               colour = "grey", alpha = 0.8, size = 0.6) +  # Gridlines for custom y-axis
  scale_x_continuous(expand = c(0, 0)) +  # Remove extra space on x-axis
  scale_y_continuous(limits = c(-50, 100)) +  # Set y-axis limits
  coord_polar() +  # Transform the plot to polar coordinates (circular)
  theme_void() +  # Remove default background and grid
  theme(legend.position = 'none')  # Remove the legend
p

# Save the plot as a PDF file
ggsave("pathway_go.pdf", p, width = 8, height = 8)
