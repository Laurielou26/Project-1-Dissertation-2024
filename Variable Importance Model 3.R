
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
plot(rf_varImp_multi, top = 20, main = "Top 20 Most Important Genes")
