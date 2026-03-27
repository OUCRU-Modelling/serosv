# Compute confidence interval for penalized_spline_model

Compute confidence interval for penalized_spline_model

## Usage

``` r
# S3 method for class 'penalized_spline_model'
compute_ci(x, ci = 0.95, ...)
```

## Arguments

- x:

  serosv models

- ci:

  confidence interval

- ...:

  arbitrary arguments

## Value

list of confidence interval for seroprevalence and foi Each confidence
interval dataframe with 4 variables, x and y for the fitted values and
ymin and ymax for the confidence interval
