# Remove spatial outliers

Removes spatial outliers using Local Moran's I statistic and moran
scatterplot.

## Usage

``` r
remove_inlier(
  x,
  y,
  ldist = 0,
  udist = 40,
  criteria = c("LM", "MP"),
  zero.policy = NULL
)
```

## Arguments

- x:

  an `sf` points object

- y:

  `character` with the name of the variable to use for depuration
  process

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
