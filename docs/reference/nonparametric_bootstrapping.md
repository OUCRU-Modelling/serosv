# Nonparametric bootstrapping for serosv model

Nonparametric bootstrapping for serosv model

## Usage

``` r
nonparametric_bootstrapping(mod, refit_func, newdat, nb = 200, ci = 0.95)
```

## Arguments

- mod:

  fitted serosv model

- refit_func:

  function to refit the model

- newdat:

  new age-range for estimating FoI

- nb:

  number of bootstrap iterations

- ci:

  confidence level for the interval
