# Load required libraries
library(cluster)
library(dplyr)
library(tibble)
library(survival)
library(survminer)
library(pheatmap)
library(ggplot2)
library(reshape2)

# Automatically install missing libraries
required_packages <- c("cluster", "dplyr", "tibble", "survival", "survminer", "pheatmap", "ggplot2", "reshape2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, library, character.only = TRUE)

# Load data (add relevant file paths)
data_clinical_patient <- read.delim("C:/Users/lauri/OneDrive - Queen's University Belfast/Master's/Dissertation/Dissertation/Dataset/Lung/lusc_tcga_pan_can_atlas_2018/data_clinical_patient.txt", comment.char="#")
data_clinical_sample <- read.delim("C:/Users/lauri/OneDrive - Queen's University Belfast/Master's/Dissertation/Dissertation/Dataset/Lung/lusc_tcga_pan_can_atlas_2018/data_clinical_sample.txt", comment.char="#")

# Data Preprocessing and Merging

# Combine clusters assignments and Clinical patient data  
sample_clusters_k_append <- sample_clusters_k

# Remove ".01" from the end of sample IDs
sample_clusters_k_append$PATIENT_ID <- sub("\\.01$", "", sample_clusters_k_append$Sample)

# Replace dots (.) with dashes (-)
sample_clusters_k_append$PATIENT_ID <- gsub("\\.", "-", sample_clusters_k_append$PATIENT_ID)

# Reset the row names of the data frame to NULL
rownames(sample_clusters_k_append) <- NULL

# Reorder the columns to make the "patient_id" column the first column
sample_clusters_k_append <- sample_clusters_k_append[, c("PATIENT_ID", names(sample_clusters_k_append)[-ncol(sample_clusters_k_append)])]

# Remove sample column
sample_clusters_k_append <- sample_clusters_k_append %>% select(-Sample_ID)

# Merge with clinical patient data
merged_data_df_K <- merge(sample_clusters_k_append, data_clinical_patient, by = "PATIENT_ID")

# Functions for statistical tests

run_chi_square_test <- function(data, var1, var2) {
  test_result <- chisq.test(table(data[[var1]], data[[var2]]))
  return(test_result$p.value)
}

run_fisher_test <- function(data, var1, var2) {
  test_result <- fisher.test(table(data[[var1]], data[[var2]]), simulate.p.value = TRUE, B = 1e5)
  return(test_result$p.value)
}

run_shapiro_wilk_test <- function(data, var) {
  # Ensure the variable is numeric
  numeric_data <- as.numeric(data[[var]])
  test_result <- shapiro.test(numeric_data)
  return(test_result$p.value)
}

run_kruskal_wallis_test <- function(data, grouping_var, var) {
  # Ensure the variable is numeric
  numeric_data <- as.numeric(data[[var]])
  test_result <- kruskal.test(numeric_data ~ data[[grouping_var]])
  return(test_result$p.value)
}

run_mann_whitney_test <- function(data, grouping_var, var) {
  numeric_data <- as.numeric(data[[var]])
  test_result <- wilcox.test(numeric_data ~ data[[grouping_var]], exact = FALSE)
  return(test_result$p.value)
}


run_tests_for_group <- function(data, grouping_var, test_vars, continuous_vars = NULL) {
  results <- lapply(test_vars, function(var) {
    p_chi <- run_chi_square_test(data, grouping_var, var)
    p_fisher <- run_fisher_test(data, grouping_var, var)
    return(c(p_chi, p_fisher))
  })
  
  if (!is.null(continuous_vars)) {
    continuous_results <- lapply(continuous_vars, function(var) {
      p_shapiro <- run_shapiro_wilk_test(data, var)
      p_kruskal <- run_kruskal_wallis_test(data, grouping_var, var)
      return(c(p_shapiro, p_kruskal))
    })
    
    results <- c(results, continuous_results)
  }
  
  return(do.call(rbind, results))
}

run_tests_for_group_cluster <- function(data, grouping_var, test_vars, continuous_vars = NULL) {
  results <- lapply(test_vars, function(var) {
    p_chi <- run_chi_square_test(data, grouping_var, var)
    p_fisher <- run_fisher_test(data, grouping_var, var)
    return(c(p_chi, p_fisher))
  })
  
  if (!is.null(continuous_vars)) {
    continuous_results <- lapply(continuous_vars, function(var) {
      p_shapiro <- run_shapiro_wilk_test(data, var)
      p_mann <- run_mann_whitney_test(data, grouping_var, var)
      return(c(p_shapiro, p_mann))
    })
    
    results <- c(results, continuous_results)
  }
  
  return(do.call(rbind, results))
}
adjust_p_values <- function(p_values) {
  p_adjusted <- p.adjust(p_values, method = "BH")
  significant <- p_adjusted < 0.05
  return(list(p_adjusted = p_adjusted, significant = significant))
}

# Data Visualization functions

# Set the desired order for clusters
merged_data_df_K$Cluster <- factor(merged_data_df_K$Cluster, levels = c("1", "2", "3", "4", "5", "6"))

# Function to generate heatmap with ordered clusters
generate_heatmap <- function(data, row_var, col_var, title) {
  heatmap_data <- table(data[[row_var]], data[[col_var]])
  
  # Reorder the clusters in the heatmap
  heatmap_data <- heatmap_data[order(as.numeric(row.names(heatmap_data))), ]
  
  pheatmap(heatmap_data, cluster_rows = FALSE, cluster_cols = FALSE, 
           display_numbers = TRUE, main = title)
}

# Generate bar plot with ordered clusters
generate_bar_plot <- function(data, x_var, y_var, fill_var, title) {
  data[[x_var]] <- factor(data[[x_var]], levels = c("1", "2", "3", "4", "5", "6"))
  
  melted_data <- melt(table(data[[x_var]], data[[y_var]]))
  
  plot_obj <- ggplot(melted_data, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = x_var, y = "Count", fill = fill_var, title = title) +
    theme_minimal()
  
  print(plot_obj)
}

# Clinical Patient Analysis

# Convert variables to factors
merged_data_df_K <- merged_data_df_K %>%
  mutate(SEX_factor = as.factor(SEX),
         Cluster = as.factor(Cluster),
         AJCC_PATHOLOGIC_TUMOR_STAGE = as.factor(AJCC_PATHOLOGIC_TUMOR_STAGE))

# Set Cluster 5 as the reference category
merged_data_df_K$Cluster <- relevel(merged_data_df_K$Cluster, ref = "5")

# Define categorical and continuous variables to test
test_vars <- c("SEX", "AJCC_PATHOLOGIC_TUMOR_STAGE", "PATH_M_STAGE", "PATH_N_STAGE", "PATH_T_STAGE", "ICD_10", "PERSON_NEOPLASM_CANCER_STATUS", "PRIOR_DX", "GENETIC_ANCESTRY_LABEL")
continuous_vars <- c("AGE")

# Run tests
test_results <- run_tests_for_group(merged_data_df_K, "Cluster", test_vars, continuous_vars)

# Adjust p-values
p_values <- unlist(test_results)
p_adjustment_results <- adjust_p_values(p_values)

# Print results
print(p_adjustment_results$p_adjusted)
print(p_adjustment_results$significant)

# Generate heatmaps and bar plots for categorical variables
for (var in test_vars) {
  generate_heatmap(merged_data_df_K, "Cluster", var, paste("Heatmap of", var, "Across Clusters"))
  generate_bar_plot(merged_data_df_K, "Cluster", var, var, paste("Distribution of", var, "Across Clusters"))
}

# Generate boxplots for continuous variables
generate_boxplot <- function(data, continuous_var, title) {
  data$Cluster <- factor(data$Cluster, levels = c("1", "2", "3", "4", "5", "6", "7"))
  
  plot_obj <- ggplot(data, aes(x = Cluster, y = get(continuous_var))) +
    geom_boxplot() +
    labs(x = "Cluster", y = continuous_var, title = title) +
    theme_minimal()
  
  print(plot_obj)
}


# Analyze clinical sample categorical variables
clinical_sample_vars <- c("TISSUE_SOURCE_SITE_CODE")


# Merge with clinical sample data
merged_data_df_K_clinical <- merge(merged_data_df_K, data_clinical_sample, by = "PATIENT_ID")

# Define continuous variables for clinical sample data
continuous_vars_clinical <- c("ANEUPLOIDY_SCORE", "MSI_SCORE_MANTIS", "TMB_NONSYNONYMOUS")

# Run tests for clinical sample variables, including continuous ones
clinical_test_results <- run_tests_for_group(merged_data_df_K_clinical, "Cluster", clinical_sample_vars, continuous_vars_clinical)
clinical_p_values <- unlist(clinical_test_results)
clinical_p_adjustment_results <- adjust_p_values(clinical_p_values)

# Print clinical sample variable results
print(clinical_p_adjustment_results$p_adjusted)
print(clinical_p_adjustment_results$significant)

merged_data_df_K_clinical$Cluster <- factor(merged_data_df_K_clinical$Cluster, levels = c("1", "2", "3", "4", "5", "6"))


# Generate heatmaps and bar plots for categorical variables in clinical sample data
for (var in clinical_sample_vars) {
  generate_heatmap(merged_data_df_K_clinical, "Cluster", var, paste("Heatmap of", var, "Across Clusters"))
  generate_bar_plot(merged_data_df_K_clinical, "Cluster", var, var, paste("Distribution of", var, "Across Clusters"))
}

# Generate boxplots for continuous variables in clinical sample data
for (var in continuous_vars_clinical) {
  plot_obj_clinical <- ggplot(merged_data_df_K_clinical, aes(x = Cluster, y = get(var))) +
    geom_boxplot() +
    labs(x = "Cluster", y = var, title = paste("Distribution of", var, "Across Clusters")) +
    theme_minimal()
  
  # Print the plot object
  print(plot_obj_clinical)
}




# Group Comparisons (Cluster 5 vs Pooled Clusters)

# For main data
group_data_df <- merged_data_df_K %>%
  mutate(Group = ifelse(Cluster == "5", "Cluster 5", "Other Clusters"))

# Run tests for groups
group_test_results <- run_tests_for_group_cluster(group_data_df, "Group", test_vars, continuous_vars)

# Adjust p-values for group comparisons
group_p_values <- unlist(group_test_results)
group_p_adjustment_results <- adjust_p_values(group_p_values)

# Print group comparison results
print(group_p_adjustment_results$p_adjusted)
print(group_p_adjustment_results$significant)

# Generate heatmaps and bar plots for group comparisons
for (var in test_vars) {
  generate_heatmap(group_data_df, "Group", var, paste("Heatmap of", var, "Across Groups"))
  generate_bar_plot(group_data_df, "Group", var, var, paste("Distribution of", var, "Across Groups"))
}

# For clinical sample data
group_data_df_clinical <- merged_data_df_K_clinical %>%
  mutate(Group = ifelse(Cluster == "5", "Cluster 5", "Other Clusters"))

# Run tests for groups
group_test_results_clinical <- run_tests_for_group_cluster(group_data_df_clinical, "Group", clinical_sample_vars, continuous_vars_clinical)

# Adjust p-values for group comparisons
group_p_values_clinical <- unlist(group_test_results_clinical)
group_p_adjustment_results_clinical <- adjust_p_values(group_p_values_clinical)

# Print group comparison results for clinical sample data
print(group_p_adjustment_results_clinical$p_adjusted)
print(group_p_adjustment_results_clinical$significant)

# Generate heatmaps and bar plots for group comparisons in clinical sample data
for (var in clinical_sample_vars) {
  generate_heatmap(group_data_df_clinical, "Group", var, paste("Heatmap of", var, "Across Groups"))
  generate_bar_plot(group_data_df_clinical, "Group", var, var, paste("Distribution of", var, "Across Groups"))
}

# Generate boxplots for continuous variables in clinical sample data
for (var in continuous_vars_clinical) {
  plot_obj_clinical_group <- ggplot(group_data_df_clinical, aes(x = Group, y = get(var))) +
    geom_boxplot() +
    labs(x = "Group", y = var, title = paste("Distribution of", var, "Across Groups")) +
    theme_minimal()
  
  # Print the plot object
  print(plot_obj_clinical_group)
}


# Generate boxplot for MSI_Score_mantis without showing outliers
plot_obj_clinical_group <- ggplot(merged_data_df_K_clinical, aes(x = Cluster, y = MSI_SCORE_MANTIS)) +
       geom_boxplot(outlier.shape = NA) +  # Hide outliers
       labs(x = "Group", y = "MSI_SCORE_MANTIS", title = "Distribution of MSI_SCORE_MANTIS Across Groups") +
      theme_minimal()+
  scale_y_continuous(limits = c(0.2, 0.4), breaks = seq(0, 0.5, by = 0.05))  # Adjust the y-axis scale, improves visualisation 

# Print the plot object
 print(plot_obj_clinical_group)
