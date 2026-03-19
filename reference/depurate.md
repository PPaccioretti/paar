# Spatial data depuration (error removal)

Filters spatial point data by removing erroneous observations based on
geometric, statistical, and spatial criteria. The function implements a
sequential depuration workflow commonly used in precision agriculture.

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

  An `sf` object with POINT geometries.

- y:

  A `character` string indicating the variable name used for filtering.
  If missing and only one attribute column is present, it is used by
  default.

- toremove:

  A `character` vector specifying which procedures to apply. Options are
  `"edges"`, `"outlier"`, and `"inlier"`. The order of execution is
  fixed and cannot be modified.

- crs:

  Coordinate reference system used when transforming longitude/latitude
  data. Can be an EPSG code or proj4string.

- buffer:

  A `numeric` value indicating the distance (in meters) for edge
  removal. Negative values are recommended to shrink boundaries.

- ylimitmax:

  Numeric upper bound for `y`. If `NA`, `Inf` is used.

- ylimitmin:

  Numeric lower bound for `y`. If `NA`, `-Inf` is used.

- sdout:

  Numeric multiplier for standard deviation used to detect global
  outliers.

- ldist:

  Numeric lower distance bound for neighborhood definition.

- udist:

  Numeric upper distance bound for neighborhood definition.

- criteria:

  Character vector specifying spatial outlier detection methods: `"LM"`
  (Local Moran) and/or `"MP"` (Moran Plot).

- zero.policy:

  Logical. If `TRUE`, allows empty neighbor sets; if `FALSE`, stops with
  an error.

- poly_border:

  Optional `sf` polygon defining field boundaries. If `NULL`, a hull is
  computed automatically.

## Value

An object of class `paar` (list) with:

- depurated_data:

  Filtered `sf` object

- condition:

  Character vector indicating the reason each observation was removed
  (or `NA` if retained)

## Details

The depuration process is applied in a fixed sequence:

1.  Edge removal (`"edges"`)

2.  Global outlier removal (`"outlier"`)

3.  Spatial outlier removal (`"inlier"`)

The `toremove` argument controls which of these steps are applied, but
\*\*does not modify the order of execution\*\*.

Available procedures are:

- edges:

  Removes points located within a specified `buffer` distance from the
  field boundary. The boundary is computed using a concave hull
  (`concaveman`) or a convex hull if the package is not available.

- outlier:

  Removes global outliers based on:

  - user-defined limits (`ylimitmin`, `ylimitmax`)

  - statistical thresholds defined as \\mean \pm sdout \times sd\\

- inlier:

  Identifies and removes spatial outliers using:

  - Local Moran's I statistic ("LM")

  - Moran scatterplot influence ("MP")

Default parameter values are tuned for precision agriculture datasets
(e.g., yield maps).

## References

Vega, A., Córdoba, M., Castro-Franco, M. et al. (2019). Protocol for
automating error removal from yield maps. *Precision Agriculture*, 20,
1030–1044.
[doi:10.1007/s11119-018-09632-8](https://doi.org/10.1007/s11119-018-09632-8)

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
