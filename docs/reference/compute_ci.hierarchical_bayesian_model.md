# Compute 95% credible interval for hierarchical Bayesian model

Return CrI for Seroprevalence and Force of Infection via parameters'
posterior distributions.

## Usage

``` r
# S3 method for class 'hierarchical_bayesian_model'
compute_ci(x, ci = 0.95, le = 100, ...)
```

## Arguments

- x:

  serosv models

- ci:

  confidence level for the interval

- ...:

  arbitrary arguments

## Value

a list of 2 data frames:

- seroprevalence estimates with columns: `x` (age), `y` (fitted
  seroprevalence), `ymin` and `ymax` (lower and upper credible interval
  bounds)

- FoI estimates with columns: `x` (age), `y` (fitted FoI), `ymin` and
  `ymax` (lower and upper credible interval bounds)
