# Remove errors from spatial data

Data can be filtered by null, edge values, global outliers and spatial
outliers or local defective observations. Default values are optimized
for precision agricultural data.

## Usage

``` r
depurate(
  x,
  y,
  toremove = c("edges", "outlier", "inlier"),
  crs = NULL,
  buffer = -10,
  ylimitmax = NA,
  ylimitmin = 0,
  sdout = 3,
  ldist = 0,
  udist = 40,
  criteria = c("LM", "MP"),
  zero.policy = NULL,
  poly_border = NULL
)
```

## Arguments

- x:

  an `sf` points object

- y:

  `character` with the name of the variable to use for
  depuration/filtering process

- toremove:

  `character` vector specifying the procedure to implement for errors
  removal. Default 'edges', 'outlier', 'inlier'. See Details.

- crs:

  coordinate reference system: integer with the EPSG code, or character
  with proj4string to convert coordinates if `x` has longitude/latitude
  data

- buffer:

  `numeric` distance in meters to be removed. Negative values are
  recommended

- ylimitmax:

  `numeric` of length 1 indicating the maximum limit for the `y`
  variable. If `NA` `Inf` is assumed

- ylimitmin:

  `numeric` of length 1 indicating the minimum limit for the `y`
  variable. If `NA` `-Inf` is assumed

- sdout:

  `numeric` values outside the interval \\mean ± sdout × sdout\\ values
  will be removed

- ldist:

  `numeric` lower distance bound to identify neighbors

- udist:

  `numeric` upper distance bound to identify neighbors

- criteria:

  `character` with "LM" and/or "MP" for methods to identify spatial
  outliers

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbors sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

- poly_border:

  `sf` object with one polygon or NULL. Can be the result of
  [`concaveman::concaveman`](https://joelgombin.github.io/concaveman/reference/concaveman.html)

## Value

an object of class `paar` with two elements:

- depurated_data:

  `sf` object with the data after the removal process

- condition:

  `character` vector with the condition of each observation

## Details

Possible values for `toremove` are one or more elements of:

- edges:

  All data points for a distance of `buffer` m from data edges are
  deleted.

- outlier:

  Values that are outside the mean±`sdout` are removed

- inlier:

  Local Moran index of spatial autocorrelation is calculated for each
  datum as a tool to identify inliers

## References

Vega, A., Córdoba, M., Castro-Franco, M. et al. Protocol for automating
error removal from yield maps. Precision Agric 20, 1030–1044 (2019).
https://doi.org/10.1007/s11119-018-09632-8

## Examples

``` r
library(sf)
data(barley, package = 'paar')
#Convert to an sf object
barley <- st_as_sf(barley, coords = c("X", "Y"), crs = 32720)

depurated <-
  depurate(barley, "Yield")
#> Concave hull algorithm is computed with
#> concavity = 2 and length_threshold = 0

# Summary of depurated data
summary(depurated)
#>       normal point             border spatial outlier MP spatial outlier LM 
#>         5673 (77%)          964 (13%)         343 (4.6%)         309 (4.2%) 
#>         global min            outlier 
#>          99 (1.3%)         6 (0.081%) 

# Keep only depurate data
depurated_data <- depurated$depurated_data
# Combine the condition for all data
all_data_condition <- cbind(depurated, barley)
```
