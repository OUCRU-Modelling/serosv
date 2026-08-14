# Compute confidence interval for Weibull model

Computes CI for Seroprevalence from model standard errors, and
(optionally) for Force of Infection via parametric bootstrap.

## Usage

``` r
# S3 method for class 'weibull_model'
compute_ci(x, ci = 0.95, le = 100, foi_ci = TRUE, ...)
```

## Arguments

- x:

  serosv models

- ci:

  confidence level for the interval

- le:

  length of age sequence for computing confidence interval, default to
  inputted age sequence for model fitting if NULL

- foi_ci:

  whether to compute CI for FoI

- ...:

  arbitrary argument

## Value

a list of 2 data frames:

- seroprevalence estimates with columns: `x` (age), `y` (fitted
  seroprevalence), `ymin` and `ymax` (lower and upper confidence
  interval bounds)

- FoI estimates with columns: `x` (age), `y` (fitted FoI), and if
  `foi_ci = TRUE`, `ymin` and `ymax` (lower and upper confidence
  interval bounds)
