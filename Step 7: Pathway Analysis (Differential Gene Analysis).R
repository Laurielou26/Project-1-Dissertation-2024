library(tidyverse)
if (!requireNamespace("ggrepel", quietly = TRUE))
  install.packages("ggrepel")
library(ggrepel)

#subset data 
filtered_rna_seq_data_50 <- read.csv("cleaned_metastatic_mrna_data.csv")

#  Prepare the data

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50[ ,-1]

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>%
  group_by(Entrez_Gene_Id) %>%
  mutate(Entrez_Gene_Id = paste(Entrez_Gene_Id, row_number(), sep = "_")) %>%
  ungroup()

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>% column_to_rownames(var = "Entrez_Gene_Id")

# Scale the data

filtered_rna_seq_data_scale <- scale(filtered_rna_seq_data_50)

# Perform k-means clustering with k = 6
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

#  'filtered_rna_seq_data_50' contains mRNA expression data
# and 'col_clusters' contains the cluster assignments for each sample

cluster_5_samples <- colnames(filtered_rna_seq_data_50)[col_clusters == 5]
other_samples <- colnames(filtered_rna_seq_data_50)[col_clusters != 5]

# Subset the data
mRNA_cluster_5 <- filtered_rna_seq_data_50[, cluster_5_samples]
mRNA_others <- filtered_rna_seq_data_50[, other_samples]

#DIFFERENTIAL EXPRESSION ANALYSIS 

if (!requireNamespace("limma", quietly = TRUE))
  install.packages("limma")
library(limma)

#  Prepare the data
mRNA_cluster_5 <- filtered_rna_seq_data_scale[, cluster_5_samples]
mRNA_others <- filtered_rna_seq_data_scale[, other_samples]

# Combine the data
all_data <- cbind(mRNA_cluster_5, mRNA_others)
condition <- factor(c(rep("Cluster_5", ncol(mRNA_cluster_5)), rep("Others", ncol(mRNA_others))))

# Create the design matrix
design <- model.matrix(~condition)
print("Design matrix:")
head(design)

#  Fit the model
fit <- lmFit(all_data, design)
fit <- eBayes(fit)

# Check the coefficient names
print("Coefficient names:")
head(colnames(fit$coefficients))

#  Get the results using the correct coefficient name
results <- topTable(fit, coef = "conditionOthers", adjust = "BH", number = Inf)

# View the first few significant DEGs
print("Top results:")
head(results)

# Filter for significant DEGs
significant_degs <- results[results$adj.P.Val < 0.05, ]

#  View results
head(significant_degs)

# Perform Pathway Enrichment Analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)
# note org.Hs.eh.Db also has a select function and this will interfere with dplyr functions. 

# List available key types
keytypes(org.Hs.eg.db)


#convert significant degs to hugo symbols 

# Remove the _1 suffix from row names
rownames(significant_degs) <- sub("_1$", "", rownames(significant_degs))

# Verify
head(rownames(significant_degs))

# Remove the _1 suffix from row names
rownames(significant_degs) <- sub("_1$", "", rownames(significant_degs))

# Verify 
head(rownames(significant_degs))

# rownames are Entrez IDs after removing suffix
entrez_ids <- rownames(significant_degs)

# Convert Entrez IDs to Gene Symbols
gene_symbols <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Merge with significant_degs to add Gene Symbols
significant_degs$ENTREZID <- rownames(significant_degs)
significant_degs_annotated <- merge(significant_degs, gene_symbols, by.x = "ENTREZID", by.y = "ENTREZID")

# Set gene symbols as rownames
rownames(significant_degs_annotated) <- significant_degs_annotated$SYMBOL

# Drop unnecessary columns
significant_degs_annotated <- significant_degs_annotated[ , !(names(significant_degs_annotated) %in% c("ENTREZID", "GENENAME"))]

# View the  results
head(significant_degs_annotated)

# Save the annotated results
write.csv(significant_degs_annotated, "significant_degs_annotated.csv")


# significant_degs$ENTREZID contains Entrez IDs of significant DEGs
entrez_ids <- rownames(significant_degs)

# Perform KEGG pathway enrichment analysis
enrich_result <- enrichKEGG(gene = entrez_ids, organism = 'hsa')

# View the results
head(enrich_result)

if (!requireNamespace("enrichplot", quietly = TRUE))
  install.packages("enrichplot")
library(enrichplot)

# Dotplot for KEGG enrichment
dotplot(enrich_result)

# Barplot for KEGG enrichment
barplot(enrich_result)


#Volcano Plots 

library(ggplot2)

# Add -log10 adjusted p-value column
significant_degs_annotated$logP <- -log10(significant_degs_annotated$adj.P.Val)

# fold change is in log2 
significant_degs_annotated$logFC <- significant_degs_annotated$logFC

# Add a column to categorize points as significant or not
significant_degs_annotated$Significant <- ifelse(significant_degs_annotated$adj.P.Val < 0.05 & abs(significant_degs_annotated$logFC) > 1, "Significant", "Not Significant")


volcano_plot <- ggplot(significant_degs_annotated, aes(x = logFC, y = logP)) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression in NSCLC Cluster 5 vs Others ", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  
  # Add labels to significant genes
  geom_text_repel(
    data = subset(significant_degs_annotated, Significant == "Significant"),
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = 20
  )

# Display the plot
print(volcano_plot)
