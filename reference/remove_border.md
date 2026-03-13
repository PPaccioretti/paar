# Remove borders

Remove borders

## Usage

``` r
remove_border(x, crs = NULL, buffer, poly_border = NULL)
```

## Arguments

- x:

  an `sf` points object

- crs:

  coordinate reference system: integer with the EPSG code, or character
  with proj4string to convert coordinates if `x` has longitude/latitude
  data

- buffer:

  `numeric` distance in meters to be removed. Negative values are
  recommended

- poly_border:

  `sf` object with one polygon or NULL. Can be the result of
  [`concaveman::concaveman`](https://joelgombin.github.io/concaveman/reference/concaveman.html)

## Details

Removes all points from `x` that are `buffer` meters from boundary.
