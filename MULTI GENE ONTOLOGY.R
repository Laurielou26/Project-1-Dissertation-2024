# Load the saved partition indices
TrainingDataIndex_multi <- readRDS("TrainingDataIndex_multi.rds")
training_samples_multi <- readRDS("training_samples_multi.rds")
test_samples_multi <- readRDS("test_samples_multi.rds")


ngenes <- ncol(training_samples_multi) - 2

dim(training_samples_multi)
dim(test_samples_multi)

set.seed(123)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 15,
  repeats = 15,
  classProbs = T,
  savePredictions = T,
  summaryFunction = multiClassSummary,
  verboseIter = TRUE,
  allowParallel = TRUE
)

# Create valid variable names for class levels
training_samples_multi$Cluster <- make.names(training_samples_multi$Cluster)

rf_multi_model2 <- train(x = training_samples_multi[1:ngenes],
                         y= training_samples_multi$Cluster,
                         method = "ranger",
                         metric = "ROC",
                         importance= 'impurity',
                         trControl = ctrl,
                         num.trees = 1000)

rf_multi_model2
rf_multi_model2 <- readRDS("rf_model_multi_cv2.rds")


#get the predictions
rfcv_pred_multi <- predict(rf_multi_model2, test_samples_multi[1:ngenes])

# Create valid variable names for class levels
test_samples_multi$Cluster <- make.names(test_samples_multi$Cluster)

# compare the predictins to the actual classifications
rfcv_cm_multi <- caret::confusionMatrix(rfcv_pred_multi, as.factor(test_samples_multi$Cluster))

rfcv_cm_multi



#Variable Importance 

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

# Plot the variable importance
plot(rf_varImp_multi, top = 20, main = "Top 20 Most Important Genes in methylation model")


#GENE ONTOLOGY 

# Install and load the required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # For human gene annotation, use the appropriate package for your species

library(clusterProfiler)
library(org.Hs.eg.db)

# List of top 20 genes
top_20_genes <- top_20_varImp_multi$Gene

# Convert gene symbols to Entrez IDs using bitr
top_20_genes_entrez <- bitr(top_20_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Identify unmapped genes
unmapped_genes <- setdiff(top_20_genes, top_20_genes_entrez$SYMBOL)

# Manually create a mapping for the unmapped genes
manual_mapping <- data.frame(
  SYMBOL = c("KRT17", "COL5A3", "POU2AF1", "GULP1", "KRT19", "TCF3", "HBEGF"), # Replace with actual gene symbols
  ENTREZID = c("3872", "50509", "5450", "51454", "3880", "6929", "1839") # Replace with corresponding Entrez IDs
)

# Combine automatic and manual mappings
complete_mapping <- rbind(top_20_genes_entrez, manual_mapping)

# Check the complete mapping
print(complete_mapping)


# Perform GO enrichment analysis
go_enrich <- enrichGO(gene         = complete_mapping$ENTREZID,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = 'ENTREZID',
                      ont          = "ALL",  # Can be "BP", "MF", "CC", or "ALL"
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

# View the results
summary(go_enrich)

# Visualize the GO enrichment results
dotplot(go_enrich, showCategory=20)

# Save the dot plot as a PNG file
ggsave("GO_enrichment_dotplot_methyl.png", plot = last_plot(), width = 10, height = 8, dpi = 300)



# Find the index of the specific gene name, e.g., "ATP2A1"
index <- which(methylation$NAME == "GULP1.1")

# Extract the corresponding row from the methylation data
methylation[index, ]


