# Compute confidence interval for the prevalence estimate from mixture model

CI for prevalence is the transform CI for mu(a) estimator, and not
accounting for the uncertainty in mu_I and mu_S estimates

## Usage

``` r
# S3 method for class 'estimate_from_mixture'
compute_ci(x, ci = 0.95, le = 100, ...)
```

## Arguments

- x:

  serosv mixture_model object

- ci:

  confidence level for the interval

- le:

  length of age sequence for computing confidence interval, default to
  inputted age sequence for model fitting if NULL

- ...:

  arbitrary arguments

## Value

a data frames of seroprevalence estimates with columns: `x` (age), `y`
(fitted seroprevalence), `ymin` and `ymax` (lower and upper confidence
interval bounds)
