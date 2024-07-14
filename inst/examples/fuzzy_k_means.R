library(sf)
data(wheat, package = 'paar')

# Transform the data.frame into a sf object
wheat_sf <- st_as_sf(wheat,
                     coords = c('x', 'y'),
                     crs = 32720)

# Run the fuzzy_k_means function
fuzzy_k_means_results <- fuzzy_k_means(wheat_sf,
                               variables = 'Tg',
                               number_cluster = 2:4)

# Print the summaryResults
fuzzy_k_means_results$summaryResults

# Print the indices
fuzzy_k_means_results$indices

# Print the cluster
head(fuzzy_k_means_results$cluster, 5)

# Combine the results in a single object
wheat_clustered <- cbind(wheat_sf, fuzzy_k_means_results$cluster)

# Plot the results
plot(wheat_clustered[, "Cluster_2"])
