# Find best smoothing parameter for Bivariate Dale model

Find best smoothing parameter for Bivariate Dale model

## Usage

``` r
find_best_bdm_param(age, y, spar_seq = seq(0, 0.1, 0.001))
```

## Arguments

- age:

  the age vector

- y:

  a dataframe of predictors (NN, NP, PN, PP), does not have to be in
  order

- spar_seq:

  sequence of smoothing parameters to be evaluated

## Value

data frame of optimal parameter for different criteria
