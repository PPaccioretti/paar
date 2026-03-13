# Summarizing paar objects

Summarizing paar objects

## Usage

``` r
# S3 method for class 'paar'
summary(object, ...)
```

## Arguments

- object:

  an object for which a summary is desired.

- ...:

  additional arguments affecting the summary produced.

## Value

An object of class summary.paar (data.frame) with the following columns:

- `condition` a character vector with the final condition.

- `n` a numeric vector with the number of rows for each condition.

- `percentage` a numeric vector with the percentage of rows for each
  condition.
