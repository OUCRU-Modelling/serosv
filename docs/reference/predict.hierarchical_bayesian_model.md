# Predict from an hierarchical bayesian model

Predict from an hierarchical bayesian model

## Usage

``` r
# S3 method for class 'hierarchical_bayesian_model'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  serosv models

- newdata:

  data.frame with age column to generate prediction

- ...:

  arbitrary arguments

## Value

list of confidence interval for seroprevalence and foi. Each confidence
interval dataframe with 4 variables, x and y for the fitted values and
ymin and ymax for the confidence interval
