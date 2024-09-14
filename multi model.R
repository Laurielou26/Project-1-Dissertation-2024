
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

#data being used


filtered_rna_seq_data_50 <- filtered_rna_seq_data_reset_50

#Remove hugo symbol from RNA data (leave entrex gene id)

filtered_rna_seq_data_50 <- filtered_rna_seq_data_50 %>% select(-Entrez_Gene_Id)

#Tranpose rna data 

rna_df <- as.data.frame(t(filtered_rna_seq_data_50[-1]))


## convert all columns from character to numeric
rna_df<- as.data.frame(apply(rna_df, 2 ,as.numeric))


## add genes ids as column names
colnames(rna_df) <- filtered_rna_seq_data_50$Hugo_Symbol


## add sample names as a column
rna_df$Sample_ID <- colnames(filtered_rna_seq_data_50[-1])


#data being used 
methylation_data<- read.delim("C:/Users/lauri/OneDrive - Queen's University Belfast/Master's/Dissertation/Dissertation/Dataset/Lung/lusc_tcga_pan_can_atlas_2018/data_methylation_hm27_hm450_merged.txt")
merged_data_df_K

methylation <- methylation_data

#remove metastatic samples

head(methylation)

# Option 1: Remove rows with NA values
methylation <- methylation[complete.cases(methylation[, -(1:4)]), ]


# Metastatic samples identified
metastatic_samples <- c("TCGA-18-3414", "TCGA-18-3417", "TCGA-33-4587", "TCGA-33-6738", 
                        "TCGA-33-A4WN", "TCGA-33-A5GW", "TCGA-33-AAS8", "TCGA-33-AASB", 
                        "TCGA-33-AASD", "TCGA-33-AASI", "TCGA-33-AASJ", "TCGA-33-AASL", 
                        "TCGA-34-8455", "TCGA-37-4132", "TCGA-43-6647", "TCGA-43-6770", 
                        "TCGA-43-6771", "TCGA-43-6773", "TCGA-43-7656", "TCGA-43-7657", 
                        "TCGA-43-8115", "TCGA-43-A56U", "TCGA-56-5897", "TCGA-56-6546", 
                        "TCGA-56-7223", "TCGA-56-7731", "TCGA-56-8082", "TCGA-56-8083", 
                        "TCGA-56-8201", "TCGA-56-8304", "TCGA-56-8308", "TCGA-56-8309", 
                        "TCGA-56-8504", "TCGA-56-8623", "TCGA-56-8624", "TCGA-56-8625", 
                        "TCGA-56-8626", "TCGA-56-8628", "TCGA-56-8629", "TCGA-56-A49D", 
                        "TCGA-56-A4BX", "TCGA-56-A4BY", "TCGA-56-A5DR", "TCGA-56-A5DS", 
                        "TCGA-56-A62T", "TCGA-58-8386", "TCGA-60-2709", "TCGA-66-2742", 
                        "TCGA-68-7756", "TCGA-68-7757", "TCGA-68-8250", "TCGA-68-A59J", 
                        "TCGA-6A-AB49", "TCGA-90-6837", "TCGA-90-7766", "TCGA-90-7767", 
                        "TCGA-90-7769", "TCGA-90-7964", "TCGA-90-A4ED", "TCGA-90-A4EE", 
                        "TCGA-90-A59Q", "TCGA-92-7340", "TCGA-92-7341", "TCGA-92-8063", 
                        "TCGA-92-8064", "TCGA-92-8065", "TCGA-94-7033", "TCGA-94-7943", 
                        "TCGA-94-8035", "TCGA-94-A5I4", "TCGA-96-7544", "TCGA-96-7545", 
                        "TCGA-J1-A4AH", "TCGA-LA-A446", "TCGA-LA-A7SW", "TCGA-MF-A522", 
                        "TCGA-NC-A5HF", "TCGA-NC-A5HP", "TCGA-NK-A5CR", "TCGA-NK-A5CX", 
                        "TCGA-O2-A52N", "TCGA-O2-A52Q", "TCGA-O2-A52S", "TCGA-O2-A52V", 
                        "TCGA-O2-A52W", "TCGA-O2-A5IB")

# Convert patient IDs to match the column names in the RNA-seq data
metastatic_columns <- paste0(gsub("-", ".", metastatic_samples), ".01")

# Check if these columns exist in the methylation dataframe
metastatic_columns <- metastatic_columns[metastatic_columns %in% colnames(methylation)]

# Remove metastatic samples from RNA-seq data
filtered_meta_methyl <- methylation %>% select(-all_of(metastatic_columns))

# Remove metastatic samples from RNA-seq data
filtered_meta_data_reset_methyl <- methylation %>% select(-all_of(metastatic_columns))

#ML model training attempt 

# Transpose the dataframe
methyl_df <- as.data.frame(t(filtered_meta_methyl))

# Make the first row as header
colnames(methyl_df) <- methyl_df[2,]

# Remove the first four row from df
methyl_df <- methyl_df[-1, ] # 4 times 
methyl_df <- methyl_df[-1, ]
methyl_df <- methyl_df[-1, ]
methyl_df <- methyl_df[-1, ]

# Reset row names to NULL
rownames(methyl_df) <- NULL

methyl_df <- t(methyl_df)

rownames_original <- rownames(methyl_df)

# Step 2: Convert the data to numeric (by converting the entire matrix)
methyl_df_numeric <- apply(methyl_df, 2, as.numeric)

# Step 3: Reassign the row names
rownames(methyl_df_numeric) <- rownames_original

# Check the structure of the new data frame
str(methyl_df_numeric)

# Convert the numeric matrix back to a data frame
methyl_df_numeric_df <- as.data.frame(methyl_df_numeric)

#add a new with median expression per gene 
methyl_df_numeric_df$med <- apply(methyl_df_numeric_df ,1 ,median)

# Calculate variance while removing NAs, 2 means function is applied to
methyl_df_numeric_df$var <- apply(methyl_df_numeric_df, 1, var)

# Filter based on variance threshold
methyl_filtered <- methyl_df_numeric_df %>% filter( med >= quantile(med,c(0.5, 0.9), na.rm =FALSE) & var >= quantile (var, c(0.5,0.9), na.rm =TRUE)) 

#remove the var and med columns 
methyl_filtered <- methyl_filtered[1:(ncol(methyl_filtered)-2)]

# View the filtered data
head(methyl_filtered)

methyl_filtered <- t(methyl_filtered)
methyl_filtered <- as.data.frame(methyl_filtered)

# Add Sample_id column with sample identifiers
methyl_filtered$Sample_ID <- colnames(filtered_meta_methyl) [-c(1, 2, 3, 4)] 


# Merge with rna_df and methyl_filtered by "Sample_ID"
merged_df <- merge(methyl_filtered, rna_df, by = "Sample_ID", suffixes = c("_methyl", "_rna"))

merged_multi_df <- merge(merged_df, sample_clusters_k, by = "Sample_ID")


# Ensure all data (except the Sample_ID column) is numeric
numeric_columns <- lapply(merged_multi_df[,-1], as.numeric)
numeric_df <- as.data.frame(numeric_columns)
numeric_df$Sample_ID <- merged_multi_df$Sample_ID

str(numeric_df$Sample_ID)

#remove NAs
numeric_df <- na.omit(numeric_df)

#make cluster a factor 

numeric_df$Cluster <- as.factor(numeric_df$Cluster)



ngenes <- ncol(numeric_df)-2

correlationMatrix <- cor(numeric_df[1:ngenes])

highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7)


#Split data into training and test sets 

#TrainingDataIndex_multi <- createDataPartition(numeric_df$Cluster, p=0.75, list = FALSE)
#training_samples_multi <- numeric_df[TrainingDataIndex_multi, -highlyCorrelated]
#test_samples_multi <- numeric_df[-TrainingDataIndex_multi, -highlyCorrelated]

# Load the saved partition indices
TrainingDataIndex_multi <- readRDS("TrainingDataIndex_multi.rds")
training_samples_multi <- readRDS("training_samples_multi.rds")
test_samples_multi <- readRDS("test_samples_multi.rds")


ngenes <- ncol(training_samples_multi) - 2

dim(training_samples_multi)
dim(test_samples_multi)

#model controls 

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

rf.model_multi_cv <- train(x = training_samples_multi[1:ngenes],
                     y= training_samples_multi$Cluster,
                     method = "ranger",
                     metric = "ROC",
                     importance= 'impurity',
                     trControl = ctrl,
                     num.trees = 1000)

rf.model_multi_cv<- readRDS("rf_model_multi_cv2.rds")
rf.model_multi_cv

#get the predictions
rfcv_pred_multi2 <- predict(rf.model_multi_cv, test_samples_multi[1:ngenes])

# Create valid variable names for class levels
test_samples_multi$Cluster <- make.names(test_samples_multi$Cluster)

# compare the predictins to the actual classifications
rfcv_cm_multi2 <- caret::confusionMatrix(rfcv_pred_multi2, as.factor(test_samples_multi$Cluster))

rfcv_cm_multi2 



#LDA MODEL 
set.seed(123)


lda.model_multi <- train(x = training_samples_multi[1:ngenes],
                   y= training_samples_multi$Cluster,
                   method = "lda",
                   metric = "ROC",
                   trControl = ctrl)

lda.model_multi <- readRDS("lda.model_multi.rds")
lda.model_multi

lda_pred <- predict(lda.model_multi, test_samples_multi[1:ngenes])
lda_pred <- factor(lda_pred, levels = levels(as.factor(test_samples_multi$Cluster)))
lda_multi_cm <- caret::confusionMatrix(lda_pred, as.factor(test_samples_multi$Cluster))
# view confusion matrix
lda_multi_cm


#KNN
set.seed(123)

knn.model_multi <- train(x = training_samples_multi[1:ngenes],
                   y= training_samples_multi$Cluster,
                   method = "kknn",
                   metric = "ROC",
                   trControl = ctrl)

knn.model_multi <- readRDS("knn_model_multi_roc.rds")
knn.model_multi

knn_pred_multi <- predict(knn.model_multi, test_samples_multi[1:ngenes])
knn_cm_multi <- caret::confusionMatrix(knn_pred_multi, as.factor(test_samples_multi$Cluster))
# view confusion matrix
knn_cm_multi


#AUC 

# Extract predictions
rf_pred <- rf.model_multi_cv$pred
lda_pred <- lda.model_multi$pred
knn_pred <- knn.model_multi$pred

# Assign model names
rf_pred$model <- "RF"
lda_pred$model <- "LDA"
knn_pred$model <- "KNN"

# Align column structures for each model's predictions
rf_pred_aligned <- rf_pred %>%
  select(pred, obs, X1:X6, rowIndex, Resample) %>%
  mutate(model = "RF")

lda_pred_aligned <- lda_pred %>%
  select(pred, obs, X1:X6, rowIndex, Resample) %>%
  mutate(model = "LDA")

knn_pred_aligned <- knn_pred %>%
  select(pred, obs, X1:X6, rowIndex, Resample) %>%
  mutate(model = "KNN")


# Combine aligned predictions
all_pred <- rbind(rf_pred_aligned, lda_pred_aligned, knn_pred_aligned)


# Plot ROC curves
g <- ggplot(all_pred, aes(m = X5, d = as.numeric(obs == "X5"), group = model, colour = model)) + 
  geom_roc(n.cuts = 0) + 
  style_roc()

print(g)

library(pROC)

# Calculate AUC for each model
rf_auc <- auc(as.numeric(rf_pred_aligned$obs == "X5"), rf_pred_aligned$X5)
lda_auc <- auc(as.numeric(lda_pred_aligned$obs == "X5"), lda_pred_aligned$X5)
knn_auc <- auc(as.numeric(knn_pred_aligned$obs == "X5"), knn_pred_aligned$X5)

aucs <- data.frame(
  model = c("RF", "LDA", "KNN"),
  AUC = c(rf_auc, lda_auc, knn_auc)
)

print(aucs)


#Variable Importance 

# Extract variable importance
rf_varImp_multi <- varImp(rf.model_multi_cv, scale = FALSE)

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


# Create a bar plot with Hugo symbols on the y-axis and importance on the x-axis
ggplot(top_20_varImp_multi, aes(x = reorder(Gene, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Flip the coordinates to have Hugo symbols on the y-axis
  labs(
    title = "Top 20 Most Important Genes",
    x = "Importance",
    y = "Hugo Symbol"
  ) +
  theme_minimal()  # Use a clean theme
