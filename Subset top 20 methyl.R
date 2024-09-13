TrainingDataIndex_m <- readRDS("TrainingDataIndex_m_save.rds")
training_samples_m <- readRDS("training_samples_m_save.rds")
test_samples_m <- readRDS("test_samples_m_save.rds")


ngenes <- ncol(training_samples_m) - 2

dim(training_samples_m)
dim(test_samples_m)

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
plot(rf_varImp_methyl, top = 20, main = "Top 20 Most Important Genes in methylation model")



# Select the top 20 genes
selected_genes <- training_samples_m[, top_20_varImp_methyl$Gene]

# Retrain the model with the selected top 20 genes
set.seed(123)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 15,
  repeats = 10,
  classProbs = T,
  savePredictions = T,
  summaryFunction = multiClassSummary,
  verboseIter = TRUE,
  allowParallel = TRUE
)

# Create valid variable names for class levels
training_samples_m$Cluster <- make.names(training_samples_m$Cluster)

rf.model_cv_top20_methyl <- train(x = selected_genes,
                           y = training_samples_m$Cluster,
                           method = "ranger",
                           metric = "ROC",
                           importance = 'impurity',
                           trControl = ctrl,
                           num.trees = 750)

# Save the model if needed
saveRDS(rf.model_cv_top20_methyl, "rf.model_cv_top20.rds")

# Predict on the test set using only the top 20 genes
rfcv_pred_top20 <- predict(rf.model_cv_top20, test_samples[, top_20_varImp$Gene])

# Compare the predictions to the actual classifications
rfcv_cm_top20 <- caret::confusionMatrix(rfcv_pred_top20, as.factor(test_samples$Cluster))

print(rfcv_cm_top20)
