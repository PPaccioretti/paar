# Removes outliers

Removes outliers

## Usage

``` r
remove_outlier(x, y, ylimitmax = NA, ylimitmin = 0, sdout = 3)
```

## Arguments

- x:

  an `sf` points object

- y:

  `character` with the name of the variable to use for depuration
  process

- ylimitmax:

  `numeric` of length 1 indicating the maximum limit for the `y`
  variable. If `NA` `Inf` is assumed

- ylimitmin:

  `numeric` of length 1 indicating the minimum limit for the `y`
  variable. If `NA` `-Inf` is assumed

- sdout:

  `numeric` values outside the interval \\mean ± sdout × sdout\\ values
  will be removed
