# Bind outlier condition to an object.

Bind outlier condition to an object.

## Usage

``` r
# S3 method for class 'paar'
cbind(..., deparse.level = 1)
```

## Arguments

- ...:

  objects to bind.

- deparse.level:

  integer controlling the construction of labels in the case of
  non-matrix-like arguments (for the default method):  
  `deparse.level = 0` constructs no labels;  
  the default `deparse.level = 1` typically and `deparse.level = 2`
  always construct labels from the argument names, see the ‘Value’
  section below.

## Value

`cbind` called with m.
