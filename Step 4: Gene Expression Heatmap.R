#load necessary libraries 
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(cluster)
library(tidyverse)
library(pheatmap)

#reset the filtered dataset 

# Ensure you have these datasets
filtered_rna_seq_data_50 <- read.csv("cleaned_metastatic_mrna_data.csv")

# Prepare the data
rna_df <- filtered_rna_seq_data_50 %>% select(-Hugo_Symbol, -Entrez_Gene_Id)

#  Scale the data and apply expression threshold filter
threshold <- 5
filtered_rna_seq_data_scale <- scale(rna_df)
highly_expressed_genes <- rowSums(abs(filtered_rna_seq_data_scale) > threshold) > 0
filtered_data <- filtered_rna_seq_data_scale[highly_expressed_genes, ]

# Log-transform the data
log_transformed_data <- log2(filtered_data + 1)

# Transpose RNA data and convert to data frame
rna_df <- as.data.frame(t(log_transformed_data))

# Add gene IDs as column names
colnames(rna_df) <- filtered_rna_seq_data_50$Entrez_Gene_Id[highly_expressed_genes]

# Ensure unique column names
colnames(rna_df) <- make.names(colnames(rna_df), unique = TRUE)

# Append cluster assignments
sample_clusters_k6 <- sample_clusters_k %>% rename(Row.names = Sample_ID)
merged_rna_df <- left_join(tibble::rownames_to_column(rna_df, var = "Row.names"), sample_clusters_k6, by = "Row.names")

# Convert row names back to sample IDs and remove the Row.names column
merged_rna_df <- column_to_rownames(merged_rna_df, var = "Row.names")

# Order the merged data frame by cluster assignments
merged_rna_df <- merged_rna_df %>% arrange(Cluster)

heatmap_data <- as.data.frame(merged_rna_df)
heatmap_data <- heatmap_data %>% arrange(Cluster)
heatmap_data <- as.data.frame(heatmap_data)
table(heatmap_data$Cluster)

# Define colors for clusters
unique_clusters <- unique(heatmap_data$Cluster)
cluster_colors <- setNames(RColorBrewer::brewer.pal(n = length(unique_clusters), name = "Set3"), unique_clusters)

# Create annotation for clusters
annotation_col <- data.frame(Cluster = factor(heatmap_data$Cluster, levels = unique_clusters))
rownames(annotation_col) <- rownames(heatmap_data)

# Define color palette for heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Generate heatmap
pheatmap(
  t(heatmap_data),
  cluster_rows = TRUE,    # Hierarchical clustering for genes
  cluster_cols = FALSE,    # 
  annotation_col = annotation_col,   # Annotation for clusters
  annotation_colors = list(Cluster = cluster_colors), # Cluster colors
  color = heatmap_colors,
  show_rownames = FALSE,  
  show_colnames = FALSE,  
  main = "RNA-Seq Heatmap (Ordered by Clusters)"
)
