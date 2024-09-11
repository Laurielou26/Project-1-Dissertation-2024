##This code will make a gene expression heatmap for the top 100 most variable genes 
### An appropriate threshold has to assigned by the user 

threshold <- 5 


# Load necessary libraries
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tibble", quietly = TRUE)) {
  install.packages("tibble")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)

#Using the appropriate filtered data from step 1
filtered_rna_seq_data_50 <- filtered_rna_seq_data_reset_50

# Step 1: Prepare the data
rna_df <- filtered_rna_seq_data_50 %>% select(-Hugo_Symbol, -Entrez_Gene_Id)

# Step 2: Scale the data and apply expression threshold filter
threshold <- 5
filtered_rna_seq_data_scale <- scale(rna_df)
highly_expressed_genes <- rowSums(abs(filtered_rna_seq_data_scale) > threshold) > 0
filtered_data <- filtered_rna_seq_data_scale[highly_expressed_genes, ]

# Step 3: Log-transform the data
log_transformed_data <- log2(filtered_data + 1)

# Step 4: Transpose RNA data and convert to data frame
rna_df <- as.data.frame(t(log_transformed_data))

# Add gene IDs as column names
colnames(rna_df) <- filtered_rna_seq_data_50$Entrez_Gene_Id[highly_expressed_genes]

# Ensure unique column names
colnames(rna_df) <- make.names(colnames(rna_df), unique = TRUE)

# Append cluster assignments
sample_clusters_k <- sample_clusters_k %>% rename(Row.names = Sample_ID)
merged_rna_df <- left_join(tibble::rownames_to_column(rna_df, var = "Row.names"), sample_clusters_k, by = "Row.names")

# Convert row names back to sample IDs and remove the Row.names column
merged_rna_df <- column_to_rownames(merged_rna_df, var = "Row.names")

# Order the merged data frame by cluster assignments
merged_rna_df <- merged_rna_df %>% arrange(Cluster)

# Select the top 100 most variable genes
top_genes <- apply(merged_rna_df %>% select(-Cluster), 1, var) %>% order(decreasing = TRUE) %>% head(100)
heatmap_data <- merged_rna_df[top_genes, ] 

#ensure Clusters are arranged in order
heatmap_data <- heatmap_data %>% arrange(Cluster)

#Transpose the data so Samples ar on the x-axis and genes are on the y-axis 
heatmap_data <- t(heatmap_data)

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
  heatmap_data,
  cluster_rows = TRUE,    # Hierarchical clustering for genes
  cluster_cols = FALSE,    # Cluster columns (samples) to improve grouping
  annotation_col = annotation_col,   # Annotation for clusters
  annotation_colors = list(Cluster = cluster_colors), # Cluster colors
  color = heatmap_colors,
  show_rownames = FALSE,  # Do not show gene names to improve visibility
  show_colnames = FALSE,  # Do not show sample names to improve visibility 
  main = "Heatmap of Top 100 Most Variable Genes"
)

# save heatmap as a PNG file with increased size
output_file <- "heatmap_gene_expression_clusters.png"
pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = list(Cluster = cluster_colors),
  color = heatmap_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Heatmap of Top 100 Most Variable Genes",
  filename = output_file,
  width = 12, height = 10
)
cat("Heatmap saved to", output_file, "\n")
