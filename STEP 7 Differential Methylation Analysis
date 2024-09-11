
methylation <- read.delim("C:/Users/lauri/OneDrive - Queen's University Belfast/Master's/Dissertation/Dissertation/Dataset/Lung/lusc_tcga_pan_can_atlas_2018/data_methylation_hm27_hm450_merged.txt")

sample_clusters_k<- read.csv("sample_clusters_k.csv")

# Assuming your original dataframe is named 'original_df'
library(tidyverse)

# Transpose the dataframe
methyl_df <- as.data.frame(t(methylation))

# Make the first row as header
colnames(methyl_df) <- methyl_df[1,]

# Remove the first row from df
methyl_df <- methyl_df[-(1:4), ]

# Reset row names to NULL
rownames(methyl_df) <- NULL

# Add Sample_id column with sample identifiers
methyl_df$Sample_ID <- colnames(methylation) [-c(1, 2, 3, 4)] # Excluding Composite.Element.REF

# Move Sample_id column to be the first column
methyl_df <- methyl_df %>% relocate(Sample_ID)


# Merge with sample_clusters_k6 by "Sample_ID"
merged_df <- merge(methyl_df, sample_clusters_k, by = "Sample_ID")

# Move Cluster column to be the first column
merged_df <- merged_df %>% relocate(Cluster)



# Save cleaned_gene_exp as a CSV file
output_file <- "merged_methylation.csv"
write.csv(merged_df, file = output_file, row.names = FALSE)

cat("Cleaned data saved to", output_file, "\n")

####################### RAN ON HPC #######################################
# Filter data to include only clusters 3 and 5
filtered_df <- merged_df %>%
  filter(Cluster %in% c(3, 5))

# Reshape data to long format
methyl_long <- filtered_df %>%
  pivot_longer(cols = -c(Sample_ID, Cluster), names_to = "Methylation", values_to = "Level")

# Ensure Cluster is a factor
methyl_long$Cluster <- as.factor(methyl_long$Cluster)


#################################RAN ON HPC#############################
#results saved in 
methyl_anova_results <- read.csv("C:/Users/lauri/OneDrive - Queen's University Belfast/Master's/Dissertation/Dissertation/Project gene exp/anova_and_shapiro_results.csv")


#tukeys FOR methyl_anova_results

# Filter rows where any significance column is TRUE
significant_df <- methyl_anova_results %>%
  filter( BH_significant == TRUE | fdr_significant == TRUE)

sig_shapiro <- significant_df %>% filter(shapiro_p_value > 0.05)

# Merge the significant results with the methylation dataset to get gene names
merged_results <- merge(significant_df, methylation, by.x = "Methylation", by.y = "ENTITY_STABLE_ID")

# Select relevant columns
final_results <- merged_results %>%
  select(Methylation, NAME, anova_p_value, BH_corrected_p, BH_significant, bonferroni_corrected_p, bonferroni_significant, fdr_corrected_p, fdr_significant)

# Print the significant entity IDs
print(final_results$NAME)

tukey_results_list <- list()

methyl_long$Cluster <- as.factor(methyl_long$Cluster)

# Loop through each significant gene and perform Tukey's HSD test
for (significant_gene in final_results$Methylation) {
  
  # Filter data for the significant gene
  subset_data <- methyl_long %>% filter(Methylation == significant_gene)
  
  # Check if there is enough data and multiple clusters to perform the test
  if (nrow(subset_data) > 0 && length(unique(subset_data$Cluster)) > 1) {
    # Perform Tukey's HSD test
    tukey_result <- TukeyHSD(aov(Level ~ Cluster, data = subset_data))
    
    # Convert Tukey's HSD result to a data frame
    tukey_df <- as.data.frame(tukey_result$Cluster)
    
    # Add a column for the gene name
    tukey_df$Methylation <- significant_gene
    
    # Store the result in the list
    tukey_results_list[[significant_gene]] <- tukey_df
  } else {
    # If not enough data or clusters, add a placeholder
    tukey_results_list[[significant_gene]] <- data.frame(
      Comparison = NA,
      Difference = NA,
      `Lower CI` = NA,
      `Upper CI` = NA,
      `P adj` = NA,
      ENTITY_STABLE_ID = significant_gene,
      gene = final_results$NAME
    )
  }
}

# Combine all results into a single data frame
tukey_results_df <- bind_rows(tukey_results_list)

# Print the combined results
#print(tukey_results_df) #COMMENTED OUT FOR APPENDING CODE 



# Merge the significant results with the methylation dataset to get gene names
merged_results <- merge(tukey_results_df, final_results, by.x = "Methylation")

# Select relevant columns
final_results <- merged_results %>%
  select(Methylation, NAME, diff, lwr, upr, "p adj", anova_p_value, BH_corrected_p, BH_significant, bonferroni_corrected_p, bonferroni_significant, fdr_corrected_p, fdr_significant)



# Sort the genes by `diff` to find hypermethylated and undermethylated genes
sorted_results <- final_results %>%
  arrange(desc(diff))

# Top 20 hypermethylated genes (highest `diff`)
top_20_hypermethylated <- sorted_results %>%
  head(20) %>%
  select(Methylation, NAME, diff, lwr, upr, `p adj`, anova_p_value, BH_corrected_p, BH_significant)

# Top 20 undermethylated genes (lowest `diff`)
top_20_undermethylated <- sorted_results %>%
  tail(20) %>%
  select(Methylation, NAME, diff, lwr, upr, `p adj`, anova_p_value, BH_corrected_p, BH_significant)

# Print the results
print("Top 20 Hypermethylated Genes:")
print(top_20_hypermethylated)

print("Top 20 Undermethylated Genes:")
print(top_20_undermethylated)



################################RAN ON HPC#################################
#REPEAT WITH CLUSTER 5 VS POOLED CLUSTERS 

# Filter data to include Cluster 5 and other clusters pooled together (1, 2, 3, 4, 6)
filtered_df_P <- merged_df %>%
  filter(Cluster %in% c(5, 1, 2, 3, 4, 6)) %>%
  mutate(Cluster = ifelse(Cluster == 5, "Cluster 5", "Other Clusters"))

# Reshape data to long format
methyl_long_P <- filtered_df_P %>%
  pivot_longer(cols = -c(Sample_ID, Cluster), names_to = "Methylation", values_to = "Level")

# Ensure Cluster is a factor
methyl_long_P$Cluster <- as.factor(methyl_long_P$Cluster)

#########################################RAN ON HPC ######################################
#RESULTS STORED IN 
methyl_anova_results_P <- read.csv("C:/Users/lauri/OneDrive - Queen's University Belfast/Master's/Dissertation/Dissertation/Project gene exp/anova_and_shapiro_results_Pooled.csv")

#tukeys FOR methyl_anova_results_P

# Filter rows where any significance column is TRUE
significant_df_p <- methyl_anova_results_P %>%
  filter( bh_significant == TRUE | fdr_significant == TRUE)

sig_shapiro_p <- significant_df_p %>% filter(shapiro_p_value > 0.05)

# Merge the significant results with the methylation dataset to get gene names
merged_results <- merge(significant_df_p, methylation, by.x = "Methylation", by.y = "ENTITY_STABLE_ID")

# Select relevant columns
final_results <- merged_results %>%
  select(Methylation, NAME, anova_p_value, bh_corrected_p, bh_significant, bonferroni_corrected_p, bonferroni_significant, fdr_corrected_p, fdr_significant)

# Print the significant entity IDs
print(final_results$NAME)

tukey_results_list <- list()

methyl_long_P$Cluster <- as.factor(methyl_long_P$Cluster)

# Loop through each significant gene and perform Tukey's HSD test
for (significant_gene in final_results$Methylation) {
  
  # Filter data for the significant gene
  subset_data <- methyl_long_P %>% filter(Methylation == significant_gene)
  
  # Check if there is enough data and multiple clusters to perform the test
  if (nrow(subset_data) > 0 && length(unique(subset_data$Cluster)) > 1) {
    # Perform Tukey's HSD test
    tukey_result <- TukeyHSD(aov(Level ~ Cluster, data = subset_data))
    
    # Convert Tukey's HSD result to a data frame
    tukey_df <- as.data.frame(tukey_result$Cluster)
    
    # Add a column for the gene name
    tukey_df$Methylation <- significant_gene
    
    # Store the result in the list
    tukey_results_list[[significant_gene]] <- tukey_df
  } else {
    # If not enough data or clusters, add a placeholder
    tukey_results_list[[significant_gene]] <- data.frame(
      Comparison = NA,
      diff = NA,
      `Lower CI` = NA,
      `Upper CI` = NA,
      `P adj` = NA,
      ENTITY_STABLE_ID = significant_gene,
      gene = final_results$NAME
    )
  }
}

# Combine all results into a single data frame
tukey_results_df <- bind_rows(tukey_results_list)

# Print the combined results
#print(tukey_results_df) #cOMMENTED OUT FOR APPENDING CODE 



# Merge the significant results with the methylation dataset to get gene names
merged_results <- merge(tukey_results_df, final_results, by.x = "Methylation")

# Select relevant columns
final_results <- merged_results %>%
  select(Methylation, NAME, diff, lwr, upr, "p adj", anova_p_value, bh_corrected_p, bh_significant, bonferroni_corrected_p, bonferroni_significant, fdr_corrected_p, fdr_significant)



# Sort the genes by `diff` to find hypermethylated and undermethylated genes
sorted_results <- final_results %>%
  arrange(desc(diff))

# Top 20 hypermethylated genes (highest `diff`)
top_20_hypermethylated_p <- sorted_results %>%
  head(20) %>%
  select(Methylation, NAME, diff, lwr, upr, `p adj`, anova_p_value, bh_corrected_p, bh_significant)

# Top 20 undermethylated genes (lowest `diff`)
top_20_undermethylated_p <- sorted_results %>%
  tail(20) %>%
  select(Methylation, NAME, diff, lwr, upr, `p adj`, anova_p_value, bh_corrected_p, bh_significant)

# Print the results
print("Top 20 Hypermethylated Genes:")
print(top_20_hypermethylated_p)

print("Top 20 Undermethylated Genes:")
print(top_20_undermethylated_p)


#LOOKING AT HYPOMETHYLATED EXPRESSION
#NOTE THE TOP 20 HYPERMETHYLATED SAMPLES ARE WITH RESPECT OF OTHER CLUSTERS TO CLUSTER 5 THEREFORE THEY ARE MORE METHYLATED IN THE OTHER CLUSTERS COMPARED TO CLUSTER 5 
# FROM THIS ANALYSIS WE FOUND THAT THERE IS A INCREASE IN EXPRESSION IN THIS GENES IN CLUSTER 5 COMPARED TO THE OTHERS, MAKES SENSE AS LESS METHYLATION = MORE EXRPRESSION 

filtered_rna_seq_data_50 <- read.csv("cleaned_metastatic_mrna_data.csv")

# Step 1: Prepare the data

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50[ ,-1]

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>%
  group_by(Entrez_Gene_Id) %>%
  mutate(Entrez_Gene_Id = paste(Entrez_Gene_Id, row_number(), sep = "_")) %>%
  ungroup()

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>% column_to_rownames(var = "Entrez_Gene_Id")

# Step 2: Scale the data

filtered_rna_seq_data_scale <- scale(filtered_rna_seq_data_50)

set.seed(123)  # For reproducibility
k6 <- 6

km_cols <- kmeans(t(filtered_rna_seq_data_scale), centers = k6, nstart = 25)
col_clusters <- km_cols$cluster

# K-means clustering on rows
km_rows <- kmeans(filtered_rna_seq_data_scale, centers = k6, nstart = 25)
row_clusters <- km_rows$cluster

# Check cluster assignments
table(row_clusters)
table(col_clusters)

# Make df for gene id and clusters
gene_clusters <- data.frame(Gene_Name = rownames(filtered_rna_seq_data_50), Cluster = row_clusters)

# Assuming `top_20_hypermethylated` contains the list of hypermethylated genes
hypomethylated_genes <- top_20_undermethylated$NAME  # Ensure this is correct

# Extract samples for Cluster 5 and Cluster 3
cluster_5_samples <- names(col_clusters[col_clusters == 5])
cluster_3_samples <- names(col_clusters[col_clusters == 3])

filtered_rna_seq_data_50 <- read.csv("cleaned_metastatic_mrna_data.csv")

# Filter mRNA data for hypermethylated genes
mRNA_data_hypo <- filtered_rna_seq_data_50[filtered_rna_seq_data_50$Hugo_Symbol %in% hypomethylated_genes, ]

# Subset mRNA data for Cluster 5 and Cluster 3
mRNA_data_cluster5 <- mRNA_data_hypo[, colnames(mRNA_data_hypo) %in% cluster_5_samples, drop = FALSE]
mRNA_data_cluster3 <- mRNA_data_hypo[, colnames(mRNA_data_hypo) %in% cluster_3_samples, drop = FALSE]

# Calculate average mRNA levels for each gene in Cluster 5 and Cluster 3
mean_mRNA_cluster5 <- rowMeans(mRNA_data_cluster5, na.rm = TRUE)
mean_mRNA_cluster3 <- rowMeans(mRNA_data_cluster3, na.rm = TRUE)

# Combine the results into a data frame for comparison
mRNA_comparison_cluster5_vs_3 <- data.frame(
  Gene = mRNA_data_hypo$Hugo_Symbol,
  Mean_mRNA_Cluster5 = mean_mRNA_cluster5,
  Mean_mRNA_Cluster3 = mean_mRNA_cluster3
)

# Add a column to show the difference
mRNA_comparison_cluster5_vs_3$Difference <- mRNA_comparison_cluster5_vs_3$Mean_mRNA_Cluster5 - mRNA_comparison_cluster5_vs_3$Mean_mRNA_Cluster3

# Print the comparison table
print(mRNA_comparison_cluster5_vs_3)

# Cluster 5 vs Others

hypomethylated_genes_p <- top_20_undermethylated_p$NAME  # Ensure this is correct

# Extract samples for Cluster 5 and all other clusters
cluster_5_samples <- names(col_clusters[col_clusters == 5])
other_clusters_samples <- names(col_clusters[col_clusters != 5])

filtered_rna_seq_data_50 <- read.csv("cleaned_metastatic_mrna_data.csv")

# Filter mRNA data for hypomethylated genes
mRNA_data_hypo_p <- filtered_rna_seq_data_50[filtered_rna_seq_data_50$Hugo_Symbol %in% hypomethylated_genes_p, ]

# Subset mRNA data for Cluster 5 and other clusters
mRNA_data_cluster5 <- mRNA_data_hypo_p[, colnames(mRNA_data_hypo_p) %in% cluster_5_samples, drop = FALSE]
mRNA_data_other_clusters <- mRNA_data_hypo_p[, colnames(mRNA_data_hypo_p) %in% other_clusters_samples, drop = FALSE]

# Calculate average mRNA levels for each gene in Cluster 5 and other clusters
mean_mRNA_cluster5 <- rowMeans(mRNA_data_cluster5, na.rm = TRUE)
mean_mRNA_other_clusters <- rowMeans(mRNA_data_other_clusters, na.rm = TRUE)

# Combine the results into a data frame for comparison
mRNA_comparison_cluster5_vs_others <- data.frame(
  Gene = mRNA_data_hypo_p$Hugo_Symbol,
  Mean_mRNA_Cluster5 = mean_mRNA_cluster5,
  Mean_mRNA_OtherClusters = mean_mRNA_other_clusters
)

# Add a column to show the difference
mRNA_comparison_cluster5_vs_others$Difference <- mRNA_comparison_cluster5_vs_others$Mean_mRNA_Cluster5 - mRNA_comparison_cluster5_vs_others$Mean_mRNA_OtherClusters

# Print the comparison table
print(mRNA_comparison_cluster5_vs_others)

#LOOK AT HYPERMETHYLATED GENES 


# Assuming `top_20_hypermethylated` contains the list of hypermethylated genes
hypermethylated_genes <- top_20_hypermethylated$NAME  # Ensure this is correct

# Extract samples for Cluster 5 and Cluster 3
cluster_5_samples <- names(col_clusters[col_clusters == 5])
cluster_3_samples <- names(col_clusters[col_clusters == 3])

filtered_rna_seq_data_50 <- read.csv("cleaned_metastatic_mrna_data.csv")

# Filter mRNA data for hypermethylated genes
mRNA_data_hyper <- filtered_rna_seq_data_50[filtered_rna_seq_data_50$Hugo_Symbol %in% hypermethylated_genes, ]

# Subset mRNA data for Cluster 5 and Cluster 3
mRNA_data_cluster5 <- mRNA_data_hyper[, colnames(mRNA_data_hyper) %in% cluster_5_samples, drop = FALSE]
mRNA_data_cluster3 <- mRNA_data_hyper[, colnames(mRNA_data_hyper) %in% cluster_3_samples, drop = FALSE]

# Calculate average mRNA levels for each gene in Cluster 5 and Cluster 3
mean_mRNA_cluster5 <- rowMeans(mRNA_data_cluster5, na.rm = TRUE)
mean_mRNA_cluster3 <- rowMeans(mRNA_data_cluster3, na.rm = TRUE)

# Combine the results into a data frame for comparison
mRNA_comparison_cluster5_vs_3 <- data.frame(
  Gene = mRNA_data_hyper$Hugo_Symbol,
  Mean_mRNA_Cluster5 = mean_mRNA_cluster5,
  Mean_mRNA_Cluster3 = mean_mRNA_cluster3
)

# Add a column to show the difference
mRNA_comparison_cluster5_vs_3$Difference <- mRNA_comparison_cluster5_vs_3$Mean_mRNA_Cluster5 - mRNA_comparison_cluster5_vs_3$Mean_mRNA_Cluster3

# Print the comparison table
print(mRNA_comparison_cluster5_vs_3)

# Cluster 5 vs Others

hypermethylated_genes_p <- top_20_hypermethylated_p$NAME  # Ensure this is correct

# Extract samples for Cluster 5 and all other clusters
cluster_5_samples <- names(col_clusters[col_clusters == 5])
other_clusters_samples <- names(col_clusters[col_clusters != 5])

filtered_rna_seq_data_50 <- read.csv("cleaned_metastatic_mrna_data.csv")

# Filter mRNA data for hypomethylated genes
mRNA_data_hyper_p <- filtered_rna_seq_data_50[filtered_rna_seq_data_50$Hugo_Symbol %in% hypermethylated_genes_p, ]

# Subset mRNA data for Cluster 5 and other clusters
mRNA_data_cluster5 <- mRNA_data_hyper_p[, colnames(mRNA_data_hyper_p) %in% cluster_5_samples, drop = FALSE]
mRNA_data_other_clusters <- mRNA_data_hyper_p[, colnames(mRNA_data_hyper_p) %in% other_clusters_samples, drop = FALSE]

# Calculate average mRNA levels for each gene in Cluster 5 and other clusters
mean_mRNA_cluster5 <- rowMeans(mRNA_data_cluster5, na.rm = TRUE)
mean_mRNA_other_clusters <- rowMeans(mRNA_data_other_clusters, na.rm = TRUE)

# Combine the results into a data frame for comparison
mRNA_comparison_cluster5_vs_others <- data.frame(
  Gene = mRNA_data_hyper_p$Hugo_Symbol,
  Mean_mRNA_Cluster5 = mean_mRNA_cluster5,
  Mean_mRNA_OtherClusters = mean_mRNA_other_clusters
)

# Add a column to show the difference
mRNA_comparison_cluster5_vs_others$Difference <- mRNA_comparison_cluster5_vs_others$Mean_mRNA_Cluster5 - mRNA_comparison_cluster5_vs_others$Mean_mRNA_OtherClusters

# Print the comparison table
print(mRNA_comparison_cluster5_vs_others)

