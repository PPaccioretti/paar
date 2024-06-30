library(sf)
data(barley, package = 'paar')
#Convert to an sf object
barley <- st_as_sf(barley,
                   coords = c("X", "Y"),
                   crs = 32720)
depurated <-
  depurate(barley,
           "Yield")

# Summary of depurated data
summary(depurated)

# Keep only depurate data
depurated_data <- depurated$depurated_data
# Combine the condition for all data
all_data_condition <- cbind(depurated, barley)
