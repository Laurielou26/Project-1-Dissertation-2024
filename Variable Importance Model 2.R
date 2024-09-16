
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
