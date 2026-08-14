# Parametric bootstrapping for serosv model

Parametric bootstrapping for serosv model

## Usage

``` r
parametric_bootstrapping(
  foi_func,
  newdat,
  coef,
  vcov,
  nb = 9999,
  alpha = 0.025
)
```

## Arguments

- foi_func:

  a function that takes coefficents, newdat, and return estimated FOI

- newdat:

  new age-range to generate FOI

- coef:

  estimated coefficients

- vcov:

  variance-covariance matrix of coefficients

- nb:

  number of bootstrap iterations

- alpha:

  significance level
