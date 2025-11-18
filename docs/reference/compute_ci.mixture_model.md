# Compute confidence interval for mixture model

Compute confidence interval for mixture model

## Usage

``` r
compute_ci.mixture_model(x, ci = 0.95, ...)
```

## Arguments

- x:

  \- serosv mixture_model object

- ci:

  \- confidence interval

- ...:

  \- arbitrary arguments

## Value

list of confidence interval for susceptible and infected. Each
confidence interval is a list with 2 items for lower and upper bound of
the interval.
