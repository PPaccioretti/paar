# Modified t test

Performs a modified t-test to assess the correlation between variables
while accounting for spatial autocorrelation. This implementation wraps
[`SpatialPack::modified.ttest`](https://rdrr.io/pkg/SpatialPack/man/modified.ttest.html).

## Usage

``` r
spatial_t_test(data, variables)
```

## Arguments

- data:

  An `sf` object containing geometry and variables, or a
  `matrix`/`data.frame` with two columns representing spatial
  coordinates (e.g., X and Y).

- variables:

  A `character` vector with the names of the variables to be tested. If
  `data` is not an `sf` object, this should be a matrix or data.frame of
  variables to test.

## Value

A `data.frame` with the following columns:

- Var1:

  Name of the first variable

- Var2:

  Name of the second variable

- corr:

  Estimated correlation coefficient

- p.value:

  P-value adjusted for spatial autocorrelation

## Details

The function computes pairwise correlations between the specified
variables and adjusts the significance test to account for spatial
dependence using coordinates. If `data` is an `sf` object, coordinates
are extracted automatically. Otherwise, coordinates must be provided as
an object with two columns.

## See also

[`modified.ttest`](https://rdrr.io/pkg/SpatialPack/man/modified.ttest.html)

## Examples

``` r
if (requireNamespace("SpatialPack", quietly = TRUE)) {
  library(sf)
  data(wheat, package = 'paar')

  # Transform the data.frame into a sf object
  wheat_sf <- st_as_sf(wheat, coords = c('x', 'y'), crs = 32720)

  # Run spatial t test
  t_test_results <-
    spatial_t_test(
      wheat_sf,
      variables = c('CE30', 'CE90')
    )

  # Print the t_test_results
  t_test_results
}
#>   Var1 Var2     corr     p.value
#> 1 CE30 CE90 0.567425 0.002973324
```
