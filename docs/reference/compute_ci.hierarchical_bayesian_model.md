# Compute 95% credible interval for hierarchical Bayesian model

Compute 95% credible interval for hierarchical Bayesian model

## Usage

``` r
# S3 method for class 'hierarchical_bayesian_model'
compute_ci(x, ...)
```

## Arguments

- x:

  serosv models

- ...:

  arbitrary arguments

## Value

list of confidence interval for seroprevalence and foi. Each confidence
interval dataframe with 4 variables, x and y for the fitted values and
ymin and ymax for the confidence interval
