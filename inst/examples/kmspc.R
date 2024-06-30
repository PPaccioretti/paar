library(sf)
data(wheat, package = 'paar')

# Transform the data.frame into a sf object
wheat_sf <- st_as_sf(wheat,
                     coords = c('x', 'y'),
                     crs = 32720)

# Run the kmspc function
kmspc_results <- kmspc(wheat_sf,
                       number_cluster = 2:4)

# Print the summaryResults
kmspc_results$summaryResults

# Print the indices
kmspc_results$indices

# Print the cluster
head(kmspc_results$cluster, 5)

# Combine the results in a single object
wheat_clustered <- cbind(wheat_sf, kmspc_results$cluster)

# Plot the results
plot(wheat_clustered[, "Cluster_2"])
