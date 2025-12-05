# Compute confidence interval for a model of serosv

Compute confidence interval for a model of serosv

## Usage

``` r
# Default S3 method
compute_ci(x, ci = 0.95, le = 100, ...)
```

## Arguments

- x:

  serosv models

- ci:

  confidence interval

- le:

  number of data for computing confidence interval

- ...:

  arbitrary argument

## Value

confidence interval dataframe with 4 variables, x and y for the fitted
values and ymin and ymax for the confidence interval
