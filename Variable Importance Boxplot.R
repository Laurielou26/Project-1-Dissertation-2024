
filtered_rna_seq_data_50 <- read.csv("cleaned_metastatic_mrna_data.csv")

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

# Remove rows with NA values
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

# ML model training attempt 

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

#  Convert the data to numeric (by converting the entire matrix)
methyl_df_numeric <- apply(methyl_df, 2, as.numeric)

# Reassign the row names
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



# Select the necessary columns: top 5 genes and the Cluster
top_genes <- c("BMS1P20", "VAMP8_methyl", "KRT5", "GULP1.1", "KRT17_rna")
plotting_data <- merged_multi_df %>% select(Sample_ID, Cluster, all_of(top_genes))

# Convert the data to long format for easier plotting with ggplot2
plotting_data_long <- melt(plotting_data, id.vars = c("Sample_ID", "Cluster"), 
                           variable.name = "Gene", value.name = "Expression")

# Load ggplot2
library(ggplot2)

# Ensure Cluster is a factor and includes all 6 levels
plotting_data_long$Cluster <- factor(plotting_data_long$Cluster, levels = 1:6)

# Generate box plots with 6 clusters
ggplot(plotting_data_long, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1) +  # Customize outliers
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +  # Facet by gene
  labs(
    title = "Expression of Top 5 Important Genes Across Clusters",
    x = "Cluster",
    y = "Expression Level"
  ) +
  theme_minimal() +  # Use a minimal theme
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5)  # Center the plot title
  )
