# Load the saved partition indices
TrainingDataIndex_multi <- readRDS("TrainingDataIndex_multi.rds")
training_samples_multi <- readRDS("training_samples_multi.rds")
test_samples_multi <- readRDS("test_samples_multi.rds")

# Extract variable importance
rf_varImp_multi <- varImp(rf_multi_model2, scale = FALSE)

# Print the variable importance
print(rf_varImp_multi)

# View the top 20 most important variables (genes)
top_20_varImp_multi <- rf_varImp_multi$importance %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  arrange(desc(Overall)) %>%
  head(20)

print(top_20_varImp_multi)


sorted_training_samples_m <- training_samples_multi[order(training_samples_multi$Cluster), ]

# Prepare heatmap data with sorted rows
heatmap_data <- as.matrix(sorted_training_samples_m[, top_20_varImp_multi$Gene])

# Prepare annotation data for the heatmap with sorted rows
annotation_data <- data.frame(Cluster = sorted_training_samples_m$Cluster)
row.names(annotation_data) <- rownames(sorted_training_samples_m)  # Ensure row names match the sample IDs in heatmap data


# Get the top 20 most important genes
top_20_varImp_with_symbols <- head(top_20_varImp_multi[order(-top_20_varImp_multi$Overall), ], 20)

# Replot with Hugo symbols as labels
# First, ensure that the `Hugo_Symbol` column is used as the row names
rownames(top_20_varImp_with_symbols) <- top_20_varImp_with_symbols$Gene


# First, make sure the order of columns in `heatmap_data` matches the order in `top_20_varImp_with_symbols`
colnames(heatmap_data) <- top_20_varImp_with_symbols$Gene

# Now plot the heatmap with Hugo symbols as column names
pheatmap(heatmap_data,
         cluster_rows = FALSE,   # Keep clusters as they are already sorted
         cluster_cols = TRUE,    # Cluster genes
         scale = "row",
         show_rownames = FALSE,  # Hide row names (samples)
         show_colnames = TRUE,   # Show column names (genes)
         annotation_row = annotation_data,  # Add cluster assignments as annotations
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  # Custom color palette
         main = "Heatmap of Top 20 Most Important Genes and Methylation Sites with Ordered Clusters",
         fontsize_col = 10,      # Adjust font size for gene names
         fontsize_row = 6)       # Adjust font size for sample labels
