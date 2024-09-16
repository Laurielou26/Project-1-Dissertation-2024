library(pheatmap)

# Load the saved partition indices
TrainingDataIndex_m <- readRDS("TrainingDataIndex_m_save.rds")
training_samples_m <- readRDS("training_samples_m_save.rds")
test_samples_m <- readRDS("test_samples_m_save.rds")


ngenes <- ncol(training_samples_m) - 2
dim(training_samples_m)
dim(test_samples_m)

#Variable Importance 

# Extract variable importance
rf_varImp_methyl <- varImp(rf.model_cv_m, scale = FALSE)

# Print the variable importance
print(rf_varImp_methyl)

# View the top 20 most important variables (genes)
top_20_varImp_methyl <- rf_varImp_methyl$importance %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  arrange(desc(Overall)) %>%
  head(20)

print(top_20_varImp_methyl)

# Plot the variable importance
plot(rf_varImp_methyl, top = 20, main = "Top 20 Most Important Genes")



sorted_training_samples_m <- training_samples_m[order(training_samples_m$Cluster), ]

# Prepare heatmap data with sorted rows
heatmap_data <- as.matrix(sorted_training_samples_m[, top_20_varImp_methyl$Gene])

# Prepare annotation data for the heatmap with sorted rows
annotation_data <- data.frame(Cluster = sorted_training_samples_m$Cluster)
row.names(annotation_data) <- rownames(sorted_training_samples_m)  # Ensure row names match the sample IDs in heatmap data


# Get the top 20 most important genes
top_20_varImp_with_symbols <- head(top_20_varImp_methyl[order(-top_20_varImp_methyl$Overall), ], 20)

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
         show_rownames = FALSE, 
         show_colnames = TRUE,   
         annotation_row = annotation_data,  # Add cluster assignments as annotations
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  
         main = "Heatmap of Top 20 Most Important Genes and Methylation Sites with Ordered Clusters",
         fontsize_col = 10,      
         fontsize_row = 6)     
