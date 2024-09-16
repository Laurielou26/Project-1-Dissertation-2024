# Ensure that 'training_samples' has 'Cluster' column with cluster assignments
library(pheatmap)
# Sort the training samples by Cluster
sorted_training_samples <- training_samples[order(training_samples$Cluster), ]

# Prepare heatmap data with sorted rows
heatmap_data <- as.matrix(sorted_training_samples[, top_20_varImp$Gene])

# Prepare annotation data for the heatmap with sorted rows
annotation_data <- data.frame(Cluster = sorted_training_samples$Cluster)
row.names(annotation_data) <- rownames(sorted_training_samples)  # Ensure row names match the sample IDs in heatmap data



#Gene names 

# Select Hugo_Symbol and Entrez_Gene_Id columns from the filtered RNA data
Hugo_symbol_df <- filtered_rna_seq_data_reset_50[, c("Hugo_Symbol", "Entrez_Gene_Id")]

# Ensure unique and valid column names for Entrez_Gene_Id
Hugo_symbol_df$Entrez_Gene_Id <- make.names(Hugo_symbol_df$Entrez_Gene_Id, unique = TRUE)

# View the resulting data frame
head(Hugo_symbol_df)

# Extract the importance data frame
varImp_df <- rf_varImp_6$importance

# Add gene names as a column in varImp_df
varImp_df$Gene <- rownames(varImp_df)

# Merge the varImp_df with Hugo_symbol_df to get the Hugo_Symbols
varImp_with_symbols <- merge(top_20_varImp, Hugo_symbol_df, by.x = "Gene", by.y = "Entrez_Gene_Id", all.x = TRUE)

# View the resulting data frame with Hugo symbols

# Get the top 20 most important genes
top_20_varImp_with_symbols <- head(varImp_with_symbols[order(-varImp_with_symbols$Overall), ], 20)

# Replot with Hugo symbols as labels
# First, ensure that the `Hugo_Symbol` column is used as the row names
rownames(top_20_varImp_with_symbols) <- top_20_varImp_with_symbols$Hugo_Symbol


# Assuming `heatmap_data` has columns as Entrez Gene IDs
# and `top_20_varImp_with_symbols` contains Entrez_Gene_Id and Hugo_Symbol columns.

# First, make sure the order of columns in `heatmap_data` matches the order in `top_20_varImp_with_symbols`
colnames(heatmap_data) <- top_20_varImp_with_symbols$Hugo_Symbol

# Now plot the heatmap with Hugo symbols as column names
pheatmap(heatmap_data,
         cluster_rows = FALSE,   # Keep clusters as they are already sorted
         cluster_cols = TRUE,    # Cluster genes
         scale = "row",
         show_rownames = FALSE, 
         show_colnames = TRUE,  
         annotation_row = annotation_data,  # Add cluster assignments as annotations
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  # Custom color palette
         main = "Heatmap of Top 20 Most Important Genes with Ordered Clusters",
         fontsize_col = 10,     
         fontsize_row = 6)       
