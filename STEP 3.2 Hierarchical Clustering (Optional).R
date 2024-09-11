
#Input optimal number of clusters from step 2 
k=6

#load necessary libraries 
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tibble", quietly = TRUE)) {
  install.packages("tibble")
}
library(cluster)
library(dplyr)
library(tibble)

#reset the filtered dataset 
filtered_rna_seq_data_50 <- filtered_rna_seq_data_reset_50

# Step 1: Prepare the data

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>% select(-Hugo_Symbol)

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>%
  group_by(Entrez_Gene_Id) %>%
  mutate(Entrez_Gene_Id = paste(Entrez_Gene_Id, row_number(), sep = "_")) %>%
  ungroup()

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>% column_to_rownames(var = "Entrez_Gene_Id")

# Step 2: Scale the data

filtered_rna_seq_data_scale <- scale(filtered_rna_seq_data_50)

set.seed(123)
hcl<-hclust(dist(t(filtered_rna_seq_data_scale)), method = "ward.D2")
clu.k6<-cutree(hcl,k)
table(clu.k6)

# Create a data frame with samples and their clusters
sample_clusters_k <- data.frame(Sample_ID = rownames(t(filtered_rna_seq_data_50)), Cluster = km$cluster)

# Output cluster sizes
table(sample_clusters_k$Cluster)

