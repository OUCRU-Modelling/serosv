# Compare models

Compare models

## Usage

``` r
compare_models(...)
```

## Arguments

- ...:

  models to be compared. Must be models created by serosv. If models'
  names are not provided, indices will be used instead for the \`model\`
  column in the returned data.frame.

## Value

a data.frame of 4 columns

- model:

  name or index of the model

- type:

  model type of the given model (a serosv model name)

- AIC:

  AIC value for the model (lower value indicates better fit)

- BIC:

  BIC value for the model (lower value indicates better fit)
