# Load the saved partition indices
TrainingDataIndex <- readRDS("TrainingDataIndex.rds")
training_samples <- readRDS("training_samples.rds")
test_samples <- readRDS("test_samples.rds")



ngenes <- ncol(training_samples) - 2

dim(training_samples)
dim(test_samples)

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
training_samples$Cluster <- make.names(training_samples$Cluster)

rf.model_cv_6<- train(x = training_samples[1:ngenes],
                      y= training_samples$Cluster,
                      method = "ranger",
                      metric = "ROC",
                      importance= 'impurity',
                      trControl = ctrl,
                      num.trees = 750)

rf.model_cv_6

rf.model_cv_6 <- readRDS("rf.model_cv_6.rds")

#predictive power 

#get the predictions
rfcv_pred_6 <- predict(rf.model_cv_6, test_samples[1:ngenes])

# Create valid variable names for class levels
test_samples$Cluster <- make.names(test_samples$Cluster)

# compare the predictins to the actual classifications
rfcv_cm_6 <- caret::confusionMatrix(rfcv_pred_6, as.factor(test_samples$Cluster))

rfcv_cm_6


#Variable Importance 

# Extract variable importance
rf_varImp_6 <- varImp(rf.model_cv_6, scale = FALSE)

# Print the variable importance
print(rf_varImp_6)

# View the top 20 most important variables (genes)
top_20_varImp <- rf_varImp_6$importance %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  arrange(desc(Overall)) %>%
  head(20)

print(top_20_varImp)

# Plot the variable importance
plot(rf_varImp_6, top = 20, main = "Top 20 Most Important Genes")


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

# Load the ggplot2 package
library(ggplot2)

# Create a bar plot with Hugo symbols on the y-axis and importance on the x-axis
ggplot(top_20_varImp_with_symbols, aes(x = reorder(Hugo_Symbol, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip the coordinates to have Hugo symbols on the y-axis
  labs(
    title = "Top 20 Most Important Genes",
    x = "Importance",
    y = "Hugo Symbol"
  ) +
  theme_minimal()  # Use a clean theme

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
training_samples_x5 <- readRDS("training_samples.rds")
test_samples_x5 <- readRDS("test_samples.rds")

# Convert the Cluster column to binary for Cluster 5 vs Others
training_samples_x5$Cluster_binary <- ifelse(training_samples_x5$Cluster == 'X5', 'Cluster_5', 'Other')
test_samples_x5$Cluster_binary <- ifelse(test_samples_x5$Cluster == 'X5', 'Cluster_5', 'Other')

# Ensure Cluster_binary is a factor
training_samples_x5$Cluster_binary <- as.factor(training_samples_x5$Cluster_binary)
test_samples_x5$Cluster_binary <- as.factor(test_samples_x5$Cluster_binary)

# Define the number of genes (features) after removing highly correlated ones
ngenes <- ncol(training_samples_x5) - 3  

# Define training control
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 4,
  classProbs = TRUE,
  savePredictions = TRUE,
  summaryFunction = twoClassSummary,  # Suitable for binary classification
  verboseIter = TRUE,
  allowParallel = TRUE
)

# Train the Random Forest model
set.seed(123)
rf_model_cv <- train(x = training_samples_x5[1:ngenes],
                     y = training_samples_x5$Cluster_binary,
                     method = "ranger",
                     metric = "ROC",  # Change to "Accuracy" if needed
                     importance = 'impurity',
                     trControl = ctrl,
                     num.trees = 500)

print(rf_model_cv)

# Predictions and evaluation
rf_pred <- predict(rf_model_cv, test_samples_x5[1:(ncol(test_samples_x5)-2)])
rf_cm <- confusionMatrix(rf_pred, as.factor(test_samples_x5$Cluster_binary))
print(rf_cm)

# Extract and plot variable importance
rf_varImp2 <- varImp(rf_model_cv, scale = FALSE)
print(rf_varImp2)

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


#Shapley Values 

library(iml)
explainer <- Predictor$new(rf_model_cv, data = test_samples_x5[1:ngenes])
shapley_values <- Shapley$new(explainer, x.interest = test_samples_x5[1, 1:ngenes])
shapley_values$plot()

