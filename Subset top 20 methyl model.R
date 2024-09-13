# Load the saved partition indices
training_samples_save <- readRDS("training_samples_m_save.rds")
test_samples_save <- readRDS("test_samples_m_save.rds")

ngenes <- ncol(training_samples_save) - 2

dim(training_samples_save)
dim(test_samples_save)



#model controls 

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
training_samples_save$Cluster <- make.names(training_samples_save$Cluster)

rf.model_m_roc <- train(x = training_samples_save[1:ngenes],
                        y= training_samples_save$Cluster,
                        method = "ranger",
                        metric = "ROC",
                        importance= 'impurity',
                        trControl = ctrl,
                        num.trees = 750)

rf.model_m_roc <- readRDS("rf_model_cv_m_roc2.rds")

#predictive power 

#get the predictions
rfcv_pred_m_roc <- predict(rf.model_m_roc, test_samples_save[1:ngenes])

# Create valid variable names for class levels
test_samples_save$Cluster <- make.names(test_samples_save$Cluster)

# compare the predictins to the actual classifications
rfcv_cm_m_roc <- caret::confusionMatrix(rfcv_pred_m_roc, as.factor(test_samples_save$Cluster))

rfcv_cm_m_roc



#Variable Importance 

# Extract variable importance
rf_varImp_methyl <- varImp(rf.model_m_roc, scale = FALSE)

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


#use ggplot for consistency 

# Load the ggplot2 package
library(ggplot2)

# Create the bar plot
ggplot(top_20_varImp_methyl, aes(x = reorder(Gene, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Top 20 Most Important Genes in Methylation Model",
    x = "Gene Name",
    y = "Importance"
  ) +
  theme_minimal() +
  coord_flip()  # Optional: flip the coordinates to make the plot horizontal


#Train model 


# Select the top 20 genes
selected_genes <- training_samples_save[, top_20_varImp_methyl$Gene]

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

training_samples_save$Cluster <- make.names(training_samples_save$Cluster)

rf.model_cv_top20_methyl <- train(x = selected_genes,
                           y = training_samples_save$Cluster,
                           method = "ranger",
                           metric = "ROC",
                           importance = 'impurity',
                           trControl = ctrl,
                           num.trees = 750)

# Save the model if needed
saveRDS(rf.model_cv_top20_methyl, "rf.model_cv_top20_methyl.rds")

test_samples_save$Cluster <- make.names(test_samples_save$Cluster)

# Predict on the test set using only the top 20 genes
rfcv_pred_top20_methyl <- predict(rf.model_cv_top20_methyl, test_samples_save[, top_20_varImp_methyl$Gene])

# Compare the predictions to the actual classifications
rfcv_cm_top20_methyl <- caret::confusionMatrix(rfcv_pred_top20_methyl, as.factor(test_samples_save$Cluster))

print(rfcv_cm_top20_methyl)
