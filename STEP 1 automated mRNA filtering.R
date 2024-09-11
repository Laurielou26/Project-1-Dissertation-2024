####USER INPUT####
#Input mrna data and clinical patient here 

data_mrna_seq_v2_rsem <- read.delim("C:/Users/lauri/OneDrive - Queen's University Belfast/Master's/Dissertation/Dissertation/Dataset/Lung/lusc_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")
data_clinical_patient <- read.delim("C:/Users/lauri/OneDrive - Queen's University Belfast/Master's/Dissertation/Dissertation/Dataset/Lung/lusc_tcga_pan_can_atlas_2018/data_clinical_patient.txt", comment.char="#") # nolint

#load necessary libraries 
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)

#Split gene expression into quartiles

mrna_data <- data_mrna_seq_v2_rsem
patient_data <- data_clinical_patient

# Remove the first column (Hugo_Symbol) and Entrez_Gene_Id column 
expression_data <- mrna_data[, -c(1, 2)]

# Calculate median expression across samples for each gene
median_expression <- apply(expression_data, 1, median)

median_expression_df <- data.frame(
  Hugo_Symbol = mrna_data$Hugo_Symbol,
  Entrez_Gene_Id = mrna_data$Entrez_Gene_Id,
  Median_Expression = median_expression
)

median_expression_sorted <- median_expression_df[order(median_expression_df$Median_Expression, decreasing = TRUE), ]

#Look at the spread of data

range_values <- range(mrna_data[, -c(1, 2)])

iqr_values <- apply(mrna_data[, -c(1, 2)], 2, IQR)
str(iqr_values)

sd_values <- apply(mrna_data[, -c(1, 2)], 2, sd)

summary_stats <- summary(mrna_data[, -c(1, 2)])


#Remove lower 50% of data (median of medians)

# Calculate the median for each gene
median_expression <- apply(expression_data, 1, median)

# Calculate the median of the medians
median_of_medians <- median(median_expression)

# Filter genes based on the median of medians threshold
filtered_data_50 <- mrna_data[median_expression > median_of_medians, ]

# Check the number of genes after filtering (should be ~half of the data)
nrow(filtered_data_50)

# Save cleaned_gene_exp as a CSV file
output_file <- "cleaned_mrna_data.csv"
write.csv(filtered_data_50, file = output_file, row.names = FALSE)

cat("Cleaned data saved to", output_file, "\n")

#remove metastatic tumour samples 


# Identify samples with metastatic tumors (PATH_M_STAGE) or Stage IV tumors (AJCC_PATHOLOGIC_TUMOR_STAGE)
metastatic_samples <- data_clinical_patient$PATIENT_ID[data_clinical_patient$PATH_M_STAGE %in% c("M1", "M1A", "M1B", "MX") | 
                                                         data_clinical_patient$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IV")]

# Print the identified metastatic samples
#print(metastatic_samples)

# Convert patient IDs to match the column names in the RNA-seq data
metastatic_columns <- paste0(gsub("-", ".", metastatic_samples), ".01")

# Check if these columns exist in the cleaned_gene_exp dataframe
metastatic_columns <- metastatic_columns[metastatic_columns %in% colnames(filtered_data_50)]

# Remove metastatic samples from RNA-seq data
filtered_rna_seq_data_50 <- dplyr::select(filtered_data_50, -all_of(metastatic_columns))

# Remove metastatic samples from RNA-seq data
filtered_rna_seq_data_reset_50 <- dplyr::select(filtered_data_50, -all_of(metastatic_columns))

# Print the first few rows of the filtered RNA-seq data to confirm the removal
head(filtered_rna_seq_data_reset_50)

# Save cleaned_gene_exp as a CSV file
output_file <- "cleaned_metastatic_mrna_data.csv"
write.csv(filtered_rna_seq_data_50, file = output_file, row.names = FALSE)

cat("Cleaned data saved to", output_file, "\n")
