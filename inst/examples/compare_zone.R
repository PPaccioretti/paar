library(sf)
data(wheat, package = "paar")

##Convert to an sf object
wheat <- sf::st_as_sf(wheat,
                      coords = c("x", "y"),
                      crs = 32720)
clusters <- paar::kmspc(
  wheat,
  variables = c('CE30', 'CE90', 'Elev', 'Pe', 'Tg'),
  number_cluster = 3:4
)
data_clusters <- cbind(wheat, clusters$cluster)
compare_zone(data_clusters,
             "Elev",
             "Cluster_3")
