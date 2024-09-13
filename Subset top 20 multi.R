# Load the saved partition indices
TrainingDataIndex_multi <- readRDS("TrainingDataIndex_multi.rds")
training_samples_multi <- readRDS("training_samples_multi.rds")
test_samples_multi <- readRDS("test_samples_multi.rds")

ngenes <- ncol(training_samples_multi) - 2

dim(training_samples_multi)
dim(test_samples_multi)


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

# Select the top 20 genes
selected_genes <- training_samples_multi[, top_20_varImp_multi$Gene]

# Retrain the model with the selected top 20 genes
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

rf.model_cv_top20_multi <- train(x = selected_genes,
                                  y = training_samples_multi$Cluster,
                                  method = "ranger",
                                  metric = "ROC",
                                  importance = 'impurity',
                                  trControl = ctrl,
                                  num.trees = 1000)

# Save the model if needed
saveRDS(rf.model_cv_top20_multi, "rf.model_cv_top20_multi.rds")

# Predict on the test set using only the top 20 genes
rfcv_pred_top20_multi <- predict(rf.model_cv_top20_multi, test_samples_multi[, top_20_varImp_multi$Gene])

# Compare the predictions to the actual classifications
rfcv_cm_top20_multi <- caret::confusionMatrix(rfcv_pred_top20_multi, as.factor(test_samples_multi$Cluster))

print(rfcv_cm_top20_multi)
