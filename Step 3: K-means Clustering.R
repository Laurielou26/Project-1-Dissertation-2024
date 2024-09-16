
#Input optimal number of clusters from Step 2 
k=6

#load necessary libraries 
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(cluster)
library(tidyverse)


#reset the filtered dataset 
filtered_rna_seq_data_50 <- read.csv("cleaned_metastatic_mrna_data.csv")

# Prepare the data

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>% select(-Hugo_Symbol)

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>%
  group_by(Entrez_Gene_Id) %>%
  mutate(Entrez_Gene_Id = paste(Entrez_Gene_Id, row_number(), sep = "_")) %>%
  ungroup()

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>% column_to_rownames(var = "Entrez_Gene_Id")

# Scale the data

filtered_rna_seq_data_scale <- scale(filtered_rna_seq_data_50)

#Perform kmeans clustering 
set.seed(123)  # For reproducibility
km <- kmeans(t(filtered_rna_seq_data_scale), centers = k, nstart = 25)

table(km$cluster)

# Create a data frame with samples and their clusters
sample_clusters_k <- data.frame(Sample_ID = rownames(t(filtered_rna_seq_data_50)), Cluster = km$cluster)

# Output cluster sizes
table(sample_clusters_k$Cluster)



# Save cleaned_gene_exp as a CSV file
output_file <- "sample_clusters_k.csv"
write.csv(sample_clusters_k, file = output_file, row.names = FALSE)

cat("Cleaned data saved to", output_file, "\n")
