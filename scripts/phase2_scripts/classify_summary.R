#install libraries
	install.packages("tidyr")
	install.packages("dplyr")
	install.packages("ggplot2")
	install.packages("readr")
#Load libraries
	library(dplyr)   # data manipulation
	library(tidyr)   # reshaping
	library(ggplot2) # visualization
	library(readr)   # load tsv file

#load and read the csv file
	mydata <- read.csv("alignment_summary.csv")

#data cleaning Inspection
	str(mydata)
	summary(mydata)
	head(mydata)

#Create Unmapped Reads
	mydata_long <- mydata %>%
  mutate(Unmapped_Reads = Total_Reads - Mapped_Reads) %>%
  pivot_longer(cols = c(Mapped_Reads, Unmapped_Reads),
               names_to = "Category",
               values_to = "Reads")

#Visualize with ggplot
	ggplot(mydata_long, aes(x = Sample, y = Reads, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Mapped vs Unmapped Reads",
       y = "Number of Reads",
       x = "Sample")

#transform the data for pie chart
	mydata_long_pct <- mydata_long %>%
  group_by(Sample) %>%
  mutate(Percent = Reads / sum(Reads) * 100)

pie chart
	ggplot(mydata_long_pct, aes(x = "", y = Percent, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~Sample) +
  theme_void() +
  labs(title = "Mapped vs Unmapped Reads (%)")



#visualizatioin for taxoomic clasification

load tsv classificaation report
	class_data = read_tsv("centrifuge_bac_class_report.tsv")

file cleaning
	str(class_data)
	summary(class_data)
	head(class_data)

# Sort the data
df_class_data <- class_data %>%
  select(name, numReads, abundance) %>%
  arrange(desc(abundance))

# Plot horizontal bar chart
ggplot(df_class_data, aes(x = reorder(name, abundance), y = abundance, fill = name)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Taxonomic Abundance",
    x = "Organism",
    y = "Abundance"
  ) +
  theme_minimal() +
  theme(legend.position = "none") 
