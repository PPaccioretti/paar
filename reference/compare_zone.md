# Compare spatial zone means

Compare spatial zone means

## Usage

``` r
compare_zone(
  data,
  variable,
  zonesCol,
  alpha = 0.05,
  join = sf::st_nearest_feature,
  returnLSD = FALSE,
  grid_dim
)
```

## Arguments

- data:

  `sf` object with zones

- variable:

  `character` or `sf` object to use for mean comparison

- zonesCol:

  `character` colname from data were zone are specified

- alpha:

  `numeric` Significance level to use for comparison

- join:

  function to use for st_join if variable is `sf` object

- returnLSD:

  `logical` when LSD calculates with spatial variance should be returned

- grid_dim:

  `numeric` grid dimentins to estimate spatial variance

## Value

`list` with differences and descriptive_stat

## References

Paccioretti, P., Córdoba, M., & Balzarini, M. (2020). FastMapping:
Software to create field maps and identify management zones in precision
agriculture. Computers and Electronics in Agriculture, 175, 105556
https://doi.org/10.1016/j.compag.2020.105556.

## Examples

``` r
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
data(wheat, package = "paar")

##Convert to an sf object
wheat <- sf::st_as_sf(wheat, coords = c("x", "y"), crs = 32720)

clusters <- paar::kmspc(
  wheat,
  variables = c('CE30', 'CE90', 'Elev', 'Pe', 'Tg'),
  number_cluster = 3:4
)

data_clusters <- cbind(wheat, clusters$cluster)

compare_zone(data_clusters, "Elev", "Cluster_3")
#> $differences
#> $differences$Elev
#>   Cluster_3     Elev                       geometry stat
#> 3         3 160.2557 MULTIPOINT ((311982.8 58006...    a
#> 1         1 160.9515 MULTIPOINT ((311962.8 58006...    b
#> 2         2 161.2599 MULTIPOINT ((312192.8 58006...    c
#> 
#> 
#> $descriptive_stat
#> Key: <n>
#>        n         Elev
#>    <int>       <char>
#> 1:  5982 160.89 (0.3)
#> 
```
