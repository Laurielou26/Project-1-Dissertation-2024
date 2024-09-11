#Step 5 Survival Analysis 

###USER INPUTS 

#select a cutoff time (months) for survival analysis
cutoff_time <- 60

#select legend labels, this needs to be the same length as the clusters (k = ?)
legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")

#Select a p-value adjustment method 
p.adjust.method = "BH"

#NOTE FOR INDIVIDUAL PAIRWISE CURVES, CLUSTER PAIRS MUST BE INPUT MANUALLY AND P-VALUES RECORDED MANUALLY. 

# Load necessary libraries
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("tibble")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
if (!requireNamespace("survminer", quietly = TRUE)) {
  install.packages("survminer")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

#Cleaned data set kmeans clustering 
library(cluster)
library(dplyr)
library(tibble)
library(survival)
library(survminer)
library(ggplot2)


#Prepare the data by combining clusters assignments and Clinical patient data  

sample_clusters_k_append <- read.csv("sample_clusters_k.csv")
data_clinical_patient <- read.delim("C:/Users/lauri/OneDrive - Queen's University Belfast/Master's/Dissertation/Dissertation/Dataset/Lung/lusc_tcga_pan_can_atlas_2018/data_clinical_patient.txt", comment.char="#") # nolint

 # sample_clusters_k variable will be needed later, so make new variable when changing

# Remove ".01" from the end of sample IDs
sample_clusters_k_append$PATIENT_ID <- sub("\\.01$", "", sample_clusters_k_append$Sample)

# Replace dots (.) with dashes (-)
sample_clusters_k_append$PATIENT_ID <- gsub("\\.", "-", sample_clusters_k_append$PATIENT_ID)


# Reset the row names of the data frame to NULL
rownames(sample_clusters_k_append) <- NULL

# Reorder the columns to make the "patient_id" column the first column
sample_clusters_k_append <- sample_clusters_k_append[, c("PATIENT_ID", names(sample_clusters_k_append)[-ncol(sample_clusters_k_append)])]

#Remove sample column
sample_clusters_k_append <- sample_clusters_k_append %>% select(-Sample_ID)

#merge with clinical patient 
merged_data_df_K <- merge(sample_clusters_k_append, data_clinical_patient, by = "PATIENT_ID")

table(merged_data_df_K$Cluster) #Ensure this is the same as sample_clusters_k output 

str(merged_data_df_K)  



#SURVIVAL ANALYSIS Overall Survival 

#create survival data frame, containing patient_id, cluster assignments, OS status and OS time
cutoff_time <- cutoff_time  
survival_data <- merged_data_df_K[merged_data_df_K$OS_MONTHS <= cutoff_time, ]

surv_OS_df <- survival_data %>%
  select(PATIENT_ID, Cluster, OS_STATUS, OS_MONTHS)

str(surv_OS_df)

#remove NAs
surv_OS_df <- na.omit(surv_OS_df)

# Convert OS_STATUS to numeric (1 = event, 0 = censored)
surv_OS_df$OS_STATUS_num <- ifelse(surv_OS_df$OS_STATUS == "0:LIVING", 0, 1)

# Create a survival object
surv_obj <- Surv(surv_OS_df$OS_MONTHS, surv_OS_df$OS_STATUS_num)

# Fit survival curves for each cluster
surv_fit <- survfit(surv_obj ~ surv_OS_df$Cluster)

# Visualize survival curves
survival_plot_all_clusters <-ggsurvplot(surv_fit, data = surv_OS_df,
                               risk.table = TRUE,
                               pval = TRUE,
                               conf.int = FALSE,
                               legend.labs = legend.labs,
                               title = "Survival Analysis by Cluster",
                               risk.table.y.text = TRUE,
                               risk.table.y.text.col = TRUE,
                               risk.table.col = "strata",
                               ggtheme = theme_bw())
survival_plot_all_clusters

# Save the plot to a PNG file
ggsave("survival_analysis_by_all_clusters.png", plot = survival_plot_all_clusters$plot, width = 10, height = 8)


## Perform pairwise log-rank tests with Bonferroni correction
pairwise_results_all_clusters <- pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS_num) ~ Cluster, data = surv_OS_df, p.adjust.method = p.adjust.method)
print(pairwise_results_all_clusters)

# Fit the Cox proportional hazards model
cox_fit_all_clusters <- coxph(Surv(OS_MONTHS, OS_STATUS_num) ~ as.factor(Cluster), data = surv_OS_df)
summary(cox_fit_all_clusters)


#INDIVIDUAL PAIRWISE SURVIVAL CURVES 
#For this code to work, selected clusters need to be changed for each combination of clusters, then the generated p-values need to be added to the original_p_values variable.

selected_clusters <- c("3", "5")
filtered_data <- surv_OS_df %>% filter(Cluster %in% selected_clusters)

# Remove rows with missing observations
filtered_data <- na.omit(filtered_data)

# Create a new survival object for the selected clusters
surv_obj_selected <- Surv(filtered_data$OS_MONTHS, filtered_data$OS_STATUS_num)

# Fit survival curves for the chosen clusters
surv_fit_selected <- survfit(surv_obj_selected ~ Cluster, data = filtered_data)

# Visualize individual pairwise survival curves for the selected clusters
individual_plots <- ggsurvplot(surv_fit_selected, data = filtered_data,
           risk.table = TRUE,
           pval = TRUE,
           conf.int = FALSE,
           legend.labs = selected_clusters,
           title = "Survival Analysis for Selected Clusters",
           risk.table.y.text = TRUE,
           risk.table.y.text.col = TRUE,
           risk.table.col = "strata",
           ggtheme = theme_bw())

individual_plots

# Save the plot to a PNG file
ggsave("survival_analysis_by_indiviual_clusters.png", plot = individual_plots$plot, width = 10, height = 8)

# Original p-values
original_p_values <- c(0.33, 0.11, 0.0011, 0.024, 0.06, 0.51, 0.38, 0.019, 0.73, 0.45, 0.73, 0.01, 0.67, 0.62, 0.33)

# Adjust p-values using the Bonferroni correction
adjusted_p_values <- p.adjust(original_p_values, method = p.adjust.method)

# Combine the adjusted p-values with the original results
adjusted_results <- data.frame(Original_P_Value = original_p_values, Adjusted_P_Value = adjusted_p_values)

# Print the adjusted results
print(adjusted_results)

# Identify significant results
pairwise_significant <- adjusted_results < 0.05

print(pairwise_significant)


#SURVIVAL 5 VS OTHER CLUSTERS POOLED TOGETHER  

# Create a new grouping variable
surv_OS_df <- surv_OS_df %>%
  mutate(Group = ifelse(Cluster == "5", "Cluster 5", "Other Clusters"))

# Create a survival object for the new grouping
surv_obj <- Surv(surv_OS_df$OS_MONTHS, surv_OS_df$OS_STATUS_num)

# Fit survival curves for the new grouping
surv_fit <- survfit(surv_obj ~ Group, data = surv_OS_df)

# Visualize survival curves
pooled_clusters <- ggsurvplot(surv_fit, data = surv_OS_df,
                    risk.table = TRUE,
                     pval = TRUE,
                     conf.int = FALSE,
                     legend.labs = c("Cluster 5", "Other Clusters"),
                     title = "Survival Analysis: Cluster 5 vs Other Clusters",
                     risk.table.y.text = TRUE,
                     risk.table.y.text.col = TRUE,
                     risk.table.col = "strata",
                     ggtheme = theme_bw())

pooled_clusters

# Save the plot to a PNG file
ggsave("survival_analysis_by_pooled_clusters.png", plot = pooled_clusters$plot, width = 10, height = 8)

# Perform log-rank test to compare the survival curves
log_rank_test <- survdiff(surv_obj ~ Group, data = surv_OS_df)

# Extract p-value from the log-rank test
p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

# Print the p-value
print(paste("P-value from log-rank test:", p_value))

# Check if the p-value is significant
significance <- p_value < 0.05
print(paste("Is the difference significant? ", significance))



#SURVIVAL ANALYSIS DISEASE FREE SURVIVAL 

surv_DFS_df <- survival_data %>%
  select(PATIENT_ID, Cluster, DFS_STATUS, DFS_MONTHS)

str(surv_DFS_df)

library(survival)
library(survminer)

surv_DFS_df <- na.omit(surv_DFS_df)

# Convert OS_STATUS to numeric (1 = event, 0 = censored)
surv_DFS_df$DFS_STATUS_num <- ifelse(surv_DFS_df$DFS_STATUS == "0:DiseaseFree", 0, 1)

# Create a survival object
surv_obj_dfs <- Surv(surv_DFS_df$DFS_MONTHS, surv_DFS_df$DFS_STATUS_num)

# Fit survival curves for each cluster
surv_fit_dfs <- survfit(surv_obj_dfs ~ surv_DFS_df$Cluster)

# Visualize survival curves
DFS_survival_curve_all<- ggsurvplot(surv_fit_dfs, data = surv_DFS_df,
                         risk.table = TRUE,
                         pval = TRUE,
                         conf.int = FALSE,
                         legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6"),
                         title = "Survival Analysis by Cluster",
                         risk.table.y.text = TRUE,
                         risk.table.y.text.col = TRUE,
                         risk.table.col = "strata",
                         ggtheme = theme_bw())

DFS_survival_curve_all

# Save the plot to a PNG file
ggsave("survival_analysis_by_all_clusters_DFS.png", plot = DFS_survival_curve_all$plot, width = 10, height = 8)

## Perform pairwise log-rank tests with BH correction
pairwise_results_dfs <- pairwise_survdiff(Surv(DFS_MONTHS, DFS_STATUS_num) ~ Cluster, data = surv_DFS_df, p.adjust.method = p.adjust.method)
print(pairwise_results_dfs)

# Fit the Cox proportional hazards model
cox_fit_dfs <- coxph(Surv(DFS_MONTHS, DFS_STATUS_num) ~ as.factor(Cluster), data = surv_DFS_df)
summary(cox_fit_dfs)

# Summarize the models
summary_OS <- summary(cox_fit_all_clusters)
summary_DFS <- summary(cox_fit_dfs)

# Extract the p-values from the summary
p_value_OS <- summary_OS$coef[,"Pr(>|z|)"]
p_value_DFS <- summary_DFS$coef[,"Pr(>|z|)"]

#SURVIVAL 5 VS OTHER CLUSTERS POOLED TOGETHER  DFS


# Create a new grouping variable
surv_DFS_df <- surv_DFS_df %>%
  mutate(Group = ifelse(Cluster == "5", "Cluster 5", "Other Clusters"))

# Create a survival object for the new grouping
surv_obj <- Surv(surv_DFS_df$DFS_MONTHS, surv_DFS_df$DFS_STATUS_num)

# Fit survival curves for the new grouping
surv_fit <- survfit(surv_obj ~ Group, data = surv_DFS_df)

# Visualize survival curves
DFS_pooled_clusters<- ggsurvplot(surv_fit, data = surv_DFS_df,
                       risk.table = TRUE,
                       pval = TRUE,
                       conf.int = FALSE,
                       legend.labs = c("Cluster 5", "Other Clusters"),
                       title = "Survival Analysis: Cluster 5 vs Other Clusters",
                       risk.table.y.text = TRUE,
                       risk.table.y.text.col = TRUE,
                       risk.table.col = "strata",
                       ggtheme = theme_bw())

DFS_pooled_clusters

# Save the plot to a PNG file
ggsave("survival_analysis_by_pooled_clusters_DFS.png", plot =DFS_pooled_clusters$plot , width = 10, height = 8)

# Perform log-rank test to compare the survival curves
log_rank_test <- survdiff(surv_obj ~ Group, data = surv_DFS_df)

# Extract p-value from the log-rank test
p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

# Print the p-value
print(paste("P-value from log-rank test:", p_value))

# Check if the p-value is significant
significance <- p_value < 0.05
print(paste("Is the difference significant? ", significance))

