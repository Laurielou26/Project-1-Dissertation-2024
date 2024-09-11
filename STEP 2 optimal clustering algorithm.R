
# specify file and parameters
#data_file <- "cleaned_metastatic_mrna_data.csv"  # Path to your filtered data file
#y <- 4  # Minimum number of clusters
#z <- 10  # Maximum number of clusters


# Load necessary libraries
if (!requireNamespace("NbClust", quietly = TRUE)) {
  install.packages("NbClust")
}
library(NbClust)

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)


# Function to get optimal clusters
get_opt_clusters <- function(data_file, y, z) {
  
  # List of indices to evaluate
  ind <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin",
           "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ball", "ptbiserial", "gap",
           "frey", "mcclain", "dunn", "sdindex")
  
  # Read data
  cleaned_gene_exp <- read.csv(data_file)
  
  # Perform data preprocessing
  cleaned_gene_exp <- cleaned_gene_exp %>%
    select(-Hugo_Symbol) %>%
    group_by(Entrez_Gene_Id) %>%
    mutate(Entrez_Gene_Id = paste(Entrez_Gene_Id, row_number(), sep = "_")) %>%
    ungroup() %>%
    column_to_rownames(var = "Entrez_Gene_Id")
  
  cleaned_gene_exp_nas <- na.omit(cleaned_gene_exp)
  
  # Initialize an empty data frame to store results
  out <- data.frame(Index = character(), Best_nc = numeric(), stringsAsFactors = FALSE)
  
  # Iterate over each index
  for (i in ind) {
    tryCatch({
      res <- NbClust(data = cleaned_gene_exp_nas, distance = "euclidean", method = "ward.D", index = i, min.nc = y, max.nc = z)
      best_nc <- res$Best.nc[1]
      b <- data.frame(Index = i, Best_nc = best_nc)
      out <- rbind(out, b)
    }, error = function(e) {
      cat("Error occurred for index", i, "Skipping...\n")
    })
  }
  
  return(out)
}


# Get optimal clusters
#results <- get_opt_clusters(data_file, y, z)

# Print results
#print(results) ##commented out for appending 

# Optionally save results to a CSV file
output_file <- "optimal_clusters_results.csv"
#write.csv(results, file = output_file, row.names = FALSE)
#cat("Results saved to", output_file, "\n")