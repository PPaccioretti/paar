# Compare means between spatial zones

Compares variable means across spatial zones using a spatially-adjusted
least significant difference (LSD) approach based on kriging variance.

The function accounts for spatial variability by estimating
semivariograms and deriving a spatial variance component, which is then
used to assess differences between zone means.

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

  an `sf` object containing the spatial zones

- variable:

  either:

  - a `character` vector with column names in `data`, or

  - an `sf` object with external variables to be compared. In this case,
    values are spatially joined to `data`.

- zonesCol:

  `character`. Column name in `data` defining zones

- alpha:

  `numeric`. Significance level for mean comparison

- join:

  function used in
  [`sf::st_join`](https://r-spatial.github.io/sf/reference/st_join.html)
  when `variable` is an external `sf` object (default:
  [`sf::st_nearest_feature`](https://r-spatial.github.io/sf/reference/st_nearest_feature.html))

- returnLSD:

  `logical`. If `TRUE`, returns the LSD value used for comparisons

- grid_dim:

  `numeric`. Grid resolution used to estimate spatial variance when
  interpolating external variables. If missing, it is automatically
  determined.

## Value

A list with:

- differences:

  list of data frames with mean comparisons per variable

- descriptive_stat:

  data frame with descriptive statistics and spatial variance

## Details

When `variable` is an external `sf` object, values are interpolated
using ordinary kriging before comparison. Otherwise, cross-validation of
the variogram model is used to estimate spatial variance.

Pairwise comparisons between zones are evaluated using a
spatially-adjusted LSD criterion:

\$\$LSD = z\_{1-\alpha/2} \times \sigma\_{spatial}\$\$

where \\\sigma\_{spatial}\\ is derived from kriging variance.

Results are presented using compact letter displays to indicate groups
of zones that are not significantly different.

## References

Paccioretti, P., Córdoba, M., & Balzarini, M. (2020). FastMapping:
Software to create field maps and identify management zones in precision
agriculture. *Computers and Electronics in Agriculture*, 175, 105556.
[doi:10.1016/j.compag.2020.105556](https://doi.org/10.1016/j.compag.2020.105556)

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
