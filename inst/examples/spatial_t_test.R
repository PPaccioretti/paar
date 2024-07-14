if (requireNamespace("SpatialPack", quietly = TRUE)) {
  library(sf)
  data(wheat, package = 'paar')

  # Transform the data.frame into a sf object
  wheat_sf <- st_as_sf(wheat,
                       coords = c('x', 'y'),
                       crs = 32720)

  # Run spatial t test
  t_test_results <-
    spatial_t_test(
      wheat_sf,
      variables = c('CE30', 'CE90'))

  # Print the t_test_results
  t_test_results


}
