# Summarize joint GWR selection

Summarize joint GWR selection

## Usage

``` r
# S3 method for class 'gwrs_selection'
summary(object, n = 10L, ...)
```

## Arguments

- object:

  A result from \[gwr_select()\].

- n:

  Number of leading candidate rows to show.

- ...:

  Additional arguments.

## Value

A summary.gwrs_selection object retaining the selected tests, local
coefficient inference/status, criterion winners, top candidates and
search status.
