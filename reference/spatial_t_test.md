# Modified t test

Performs a modified version of the t test to assess the correlation
between spatial processes. See SpatialPack::modified.ttest for details.

## Usage

``` r
spatial_t_test(data, variables)
```

## Arguments

- data:

  `sf` data to extract coordinates or two columns `matrix` or
  `data.frame` specifying coordinates.

- variables:

  `character` vector with column names to perform ttest

## Value

a data.frame with the correlation and p-value for each pair of variables
