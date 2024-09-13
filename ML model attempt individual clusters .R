#Load libraries 

library(tidyverse)
library(caret) 
library(Rcpp)
library(plotly)
library(ggfortify)
library(kknn)
library(reshape2)
library(Rcpp)
library(doParallel)
library(survival)
library(survminer)
library(plotROC)

#ML model training attempt 

#data being used

filtered_rna_seq_data_50 <- filtered_rna_seq_data_reset_50
merged_data_df_K

#visualise data 

table(merged_data_df_K$PATIENT_ID, merged_data_df_K$Cluster_binary)

#Remove hugo symbol from RNA data (leave entrex gene id)

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>% select(-Hugo_Symbol)

#I already filtered the data previously, now remove genes with low variance? This removes alot of data, no further filtering will be applied 

#filtered_rna_seq_data_50$var <- apply(filtered_rna_seq_data_50[-1], 1, var)

#rna_filtered <- filtered_rna_seq_data_50 %>% filter( var >= quantile(var, c(0.5, 0.9))[2])

#rna_filtered <- rna_filtered[1:(ncol(rna_filtered)-2)]


#Tranpose rna data 

rna_df <- as.data.frame(t(filtered_rna_seq_data_50[-1]))


## convert all columns from character to numeric
rna_df<- as.data.frame(apply(rna_df, 2 ,as.numeric))


## add genes ids as column names
colnames(rna_df) <- filtered_rna_seq_data_50$Entrez_Gene_Id


## add sample names as a column
rna_df$sample <- colnames(filtered_rna_seq_data_50[-1])

# Move Sample column to be the first column
#rna_df <- rna_df %>% relocate(sample)

# Ensure unique column names
colnames(rna_df) <- make.names(colnames(rna_df), unique = TRUE)

## append cluster assignments  
merged_rna_df <- left_join(rna_df, sample_clusters_k6, by = c("sample" = "Sample_ID"))

# View the merged data frame
head(merged_rna_df[1:5]) 

# Move cluster column to be the first column
#merged_rna_df <- merged_rna_df %>% relocate(Cluster)

# Convert cluster label to binary (1 for cluster 5, 0 for others)
#merged_rna_df$Cluster_binary <- ifelse(merged_rna_df$Cluster == 5, "Cluster 5", "Other Clusters")

# Convert Cluster_binary to a factor
#merged_rna_df$Cluster_binary <- factor(merged_rna_df$Cluster_binary, levels = c(0, 1))

# Move cluster binary column to be the first column
#merged_rna_df <- merged_rna_df %>% relocate(Cluster_binary)

#make cluster a factor 

merged_rna_df$Cluster <- as.factor(merged_rna_df$Cluster)


#now need to remove highly correlated genes 


ngenes <- ncol(merged_rna_df)-2

correlationMatrix <- cor(merged_rna_df[1:ngenes])

highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7)


#Split data into training and test sets 

TrainingDataIndex <- createDataPartition(merged_rna_df$Cluster, p=0.75, list = FALSE)
training_samples <- merged_rna_df[TrainingDataIndex, -highlyCorrelated]
test_samples <- merged_rna_df[-TrainingDataIndex, -highlyCorrelated]



ngenes <- ncol(training_samples) - 2

dim(training_samples)
dim(test_samples)

#model controls 

set.seed(123)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  classProbs = T,
  savePredictions = T,
  summaryFunction = multiClassSummary,
  verboseIter = TRUE,
  allowParallel = TRUE
)

# Create valid variable names for class levels
training_samples$Cluster <- make.names(training_samples$Cluster)

rf.model_cv <- train(x = training_samples[1:ngenes],
                     y= training_samples$Cluster,
                     method = "ranger",
                     metric = "ROC",
                     importance= 'impurity',
                     trControl = ctrl,
                     num.trees = 500)

rf.model_cv


#prdictive power 

#get the predictions
rfcv_pred <- predict(rf.model_cv, test_samples[1:ngenes])

# Create valid variable names for class levels
test_samples$Cluster <- make.names(test_samples$Cluster)

# compare the predictins to the actual classifications
rfcv_cm <- caret::confusionMatrix(rfcv_pred, as.factor(test_samples$Cluster))

rfcv_cm



#LDA MODEL 
set.seed(123)


lda.model <- train(x = training_samples[1:ngenes],
                   y= training_samples$Cluster,
                   method = "lda",
                   metric = "ROC",
                   trControl = ctrl)

lda_pred <- predict(lda.model, test_samples[1:ngenes])
lda_pred <- factor(lda_pred, levels = levels(as.factor(test_samples$Cluster)))
lda_cm <- caret::confusionMatrix(lda_pred, as.factor(test_samples$Cluster))
# view confusion matrix
lda_cm


#KNN
set.seed(123)

knn.model <- train(x = training_samples[1:ngenes],
                   y= training_samples$Cluster,
                   method = "kknn",
                   metric = "ROC",
                   trControl = ctrl)


knn_pred <- predict(knn.model, test_samples[1:ngenes])
knn_cm <- caret::confusionMatrix(knn_pred, as.factor(test_samples$Cluster))
# view confusion matrix
knn_cm


#COMPARE MODELS 
# Combine predictions from all models
rf_roc <- rf.model_cv$pred[1:8]
rf_roc$model <- "RF"

lda_roc <- lda.model$pred[1:8]
lda_roc$model <- "LDA"

knn_roc <- knn.model$pred[1:8]
knn_roc$model <- "KNN"

# Combine the aligned data frames
all_roc <- rbind(rf_roc, lda_roc, knn_roc)


g <- ggplot(all_roc, aes(m="X5", d=factor(obs, levels = c("X1", "X2", "X3", "X4", "X5", "X6")),group=model, colour=model)) + 
  geom_roc(n.cuts=0) + 
  style_roc()

g

# calculate area under curve
calc_auc(g)


#Gradiant boosting 

# Gradient Boosting Machines (GBM)
gbmGrid <- expand.grid(interaction.depth = c(1, 5, 9),
                       n.trees = c(100, 500, 1000),
                       shrinkage = c(0.01, 0.1),
                       n.minobsinnode = 10)

gbm.model <- train(x = training_samples[1:ngenes],
                   y = training_samples$Cluster,
                   method = "gbm",
                   metric = "ROC",
                   tuneGrid = gbmGrid,
                   trControl = ctrl,
                   verbose = FALSE)

# Predictive power
gbm_pred <- predict(gbm.model, test_samples[1:ngenes])
gbm_cm <- caret::confusionMatrix(gbm_pred, as.factor(test_samples$Cluster))
gbm_cm



# Support Vector Machines (SVM)
svmGrid <- expand.grid(C = c(0.1, 1, 10),
                       sigma = c(0.01, 0.05, 0.1))

svm.model <- train(x = training_samples[1:ngenes],
                   y = training_samples$Cluster,
                   method = "svmRadial",
                   metric = "ROC",
                   tuneGrid = svmGrid,
                   trControl = ctrl)

# Predictive power
svm_pred <- predict(svm.model, test_samples[1:ngenes])
svm_cm <- caret::confusionMatrix(svm_pred, as.factor(test_samples$Cluster))
svm_cm

# Extreme Gradient Boosting (XGBoost)
xgbGrid <- expand.grid(nrounds = c(100, 200),
                       max_depth = c(3, 6, 9),
                       eta = c(0.01, 0.1, 0.3),
                       gamma = 0,
                       colsample_bytree = 1,
                       min_child_weight = 1,
                       subsample = 1)

xgb.model <- train(x = training_samples[1:ngenes],
                   y = training_samples$Cluster,
                   method = "xgbTree",
                   metric = "ROC",
                   tuneGrid = xgbGrid,
                   trControl = ctrl)

# Predictive power
xgb_pred <- predict(xgb.model, test_samples[1:ngenes])
xgb_cm <- caret::confusionMatrix(xgb_pred, as.factor(test_samples$Cluster))
xgb_cm


# Neural Networks (NN)
nnGrid <- expand.grid(size = c(5, 10), decay = c(0.1, 0.5))

nn.model <- train(
  x = training_samples[1:ngenes],
  y = training_samples$Cluster,
  method = "nnet",
  metric = "ROC",
  tuneGrid = nnGrid,
  trControl = ctrl,
  linout = FALSE
)


# Predictive power
nn_pred <- predict(nn.model, test_samples[1:ngenes])
nn_cm <- caret::confusionMatrix(nn_pred, as.factor(test_samples$Cluster))
nn_cm


# Combine ROC data from all models
gbm_roc <- gbm.model$pred[1:4]
gbm_roc$model <- "GBM"
svm_roc <- svm.model$pred[1:4]
svm_roc$model <- "SVM"
xgb_roc <- xgb.model$pred[1:4]
xgb_roc$model <- "XGBoost"


all_roc_1 <- rbind(all_roc, gbm_roc, svm_roc, xgb_roc)

# Plot the ROC curves
g <- ggplot(all_roc_1, aes(m = Cluster.5, d = factor(obs, levels = c("Cluster.5", "Other.Clusters")), group = model, colour = model)) + 
  geom_roc(n.cuts = 0) + 
  style_roc()

g

# Calculate Area Under Curve
calc_auc(g)

#Try to improve accuracy 

install.packages("caret")
install.packages("xgboost")
install.packages("remotes")
remotes::install_github("cran/DMwR")


# For SMOTE

library(caret)
library(xgboost)
library(DMwR)  # Provides SMOTE function

# Assuming training_samples and test_samples are already defined

# Convert Cluster_binary to factor
training_samples$Cluster_binary <- as.factor(training_samples$Cluster_binary)
test_samples$Cluster_binary <- as.factor(test_samples$Cluster_binary)

ctrl <- trainControl(
  method = "cv",
  number = 10,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  verboseIter = TRUE,
  sampling = "smote"  # Apply SMOTE to balance classes
)

xgbGrid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 6, 9),
  eta = c(0.01, 0.1, 0.3),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

xgb.model <- train(
  x = training_samples[1:ngenes],
  y = training_samples$Cluster_binary,
  method = "xgbTree",
  metric = "ROC",
  tuneGrid = xgbGrid,
  trControl = ctrl
)

# Predict on test data
xgb_pred <- predict(xgb.model, test_samples[1:ngenes])

# Generate confusion matrix
xgb_cm <- caret::confusionMatrix(xgb_pred, test_samples$Cluster_binary)
print(xgb_cm)

