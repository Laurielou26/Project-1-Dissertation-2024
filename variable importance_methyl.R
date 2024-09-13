# Load the saved partition indices
TrainingDataIndex_m <- readRDS("TrainingDataIndex_m.rds")
training_samples_m <- readRDS("training_samples_m.rds")
test_samples_m <- readRDS("test_samples_m.rds")


ngenes <- ncol(training_samples_m) - 2

dim(training_samples_m)
dim(test_samples_m)


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

rf.model_cv_m1 <- train(x = training_samples_m[1:ngenes],
                        y= training_samples_m$Cluster,
                        method = "ranger",
                        metric = "ROC",
                        importance= 'impurity',
                        trControl = ctrl,
                        num.trees = 750)

rf.model_cv_m1

rf.model_cv_m1 <- readRDS("rf_model_cv_m1.rds")

rf.model_cv_m1



#get the predictions
rfcv_pred_m1 <- predict(rf.model_cv_m1, test_samples_m[1:ngenes])

# Create valid variable names for class levels
test_samples_m$Cluster <- make.names(test_samples_m$Cluster)

# compare the predictins to the actual classifications
rfcv_cm_m1 <- caret::confusionMatrix(rfcv_pred_m1, as.factor(test_samples_m$Cluster))

rfcv_cm_m1





#Variable Importance 

# Extract variable importance
rf_varImp_methyl <- varImp(rf.model_cv_m1, scale = FALSE)

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


#Partial Dependence Plot 
library(pdp)
partial_plot <- partial(rf.model_cv_6, pred.var = "X3852", plot = TRUE)

#Using vivid package to visualise variable importance 

if (!requireNamespace("graph", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("graph")
}
install.packages("zenplots")
install.packages("vivid")

library(vivid)
library(randomForest) # if using randomForest
library(ranger)       # if using ranger

rf.model_cv_6

# Ensure you have the response variable and training data
response_variable <- "Cluster" # Adjust if necessary
training_data <- training_samples

# Calculate variable importance and interactions
viviRf <- vivi(
  fit = rf.model_cv_6$finalModel,  # Extract the final model from caret
  data = training_data,
  response = response_variable,
  gridSize = 50,
  importanceType = "agnostic",
  nmax = 500
)

viviHeatmap(mat = viviRf)

#Class specific importance X5

# Load libraries
library(caret)
library(ranger)
library(dplyr)

# Load the training and test data frames
training_samples_x5_methyl <- readRDS("training_samples_m.rds")
test_samples_x5_methyl <- readRDS("test_samples_m.rds")

# Convert the Cluster column to binary for Cluster 5 vs Others
training_samples_x5_methyl$Cluster_binary <- ifelse(training_samples_x5_methyl$Cluster == 'X5', 'Cluster_5', 'Other')
test_samples_x5_methyl$Cluster_binary <- ifelse(test_samples_x5_methyl$Cluster == 'X5', 'Cluster_5', 'Other')

# Ensure Cluster_binary is a factor
training_samples_x5_methyl$Cluster_binary <- as.factor(training_samples_x5_methyl$Cluster_binary)
test_samples_x5_methyl$Cluster_binary <- as.factor(test_samples_x5_methyl$Cluster_binary)

# Define the number of genes (features) after removing highly correlated ones
ngenes <- ncol(training_samples_x5_methyl) - 3  

# Define training control
ctrl <- trainControl(
  method = "repeatedcv",
  number = 15,
  repeats = 10,
  classProbs = TRUE,
  savePredictions = TRUE,
  summaryFunction = twoClassSummary,  # Suitable for binary classification
  verboseIter = TRUE,
  allowParallel = TRUE
)

# Train the Random Forest model
set.seed(123)
rf_model_cv_methyl <- train(x = training_samples_x5_methyl[1:ngenes],
                     y = training_samples_x5_methyl$Cluster_binary,
                     method = "ranger",
                     metric = "ROC",  # Change to "Accuracy" if needed
                     importance = 'impurity',
                     trControl = ctrl,
                     num.trees = 750)

print(rf_model_cv_methyl)

# Predictions and evaluation
rf_pred_methyl <- predict(rf_model_cv_methyl, test_samples_x5_methyl[1:(ncol(test_samples_x5_methyl)-2)])
rf_cm_methyl <- confusionMatrix(rf_pred_methyl, as.factor(test_samples_x5_methyl$Cluster_binary))
print(rf_cm_methyl)

# Extract and plot variable importance
rf_varImp2_methyl <- varImp(rf_model_cv_methyl, scale = FALSE)
print(rf_varImp2_methyl)

# View and plot top important genes
top_20_varImp2 <- rf_varImp2$importance %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  arrange(desc(Overall)) %>%
  head(20)

print(top_20_varImp2)

# Merge with gene symbols for plotting
Hugo_symbol_df <- filtered_rna_seq_data_50 %>% select(Hugo_Symbol, Entrez_Gene_Id)

# Ensure unique and valid column names for Entrez_Gene_Id
Hugo_symbol_df$Entrez_Gene_Id <- make.names(Hugo_symbol_df$Entrez_Gene_Id, unique = TRUE)

# View the resulting data frame
head(Hugo_symbol_df)
varImp_with_symbols2 <- merge(top_20_varImp2, Hugo_symbol_df, by.x = "Gene", by.y = "Entrez_Gene_Id", all.x = TRUE)

# Plot
library(ggplot2)
ggplot(varImp_with_symbols2, aes(x = reorder(Hugo_Symbol, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(
    title = "Top 20 Most Important Genes for Class X5",
    x = "Importance",
    y = "Hugo Symbol"
  ) +
  theme_minimal()

#Partial Dependence Plot 
library(pdp)
partial_plot <- partial(rf_model_cv, pred.var = "Gene_of_interest", plot = TRUE)

