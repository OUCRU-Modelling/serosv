# Compare models

Compare models

## Usage

``` r
compare_models(data, method = "AIC/BIC", ...)
```

## Arguments

- data:

  input data to fit into the models

- method:

  method to compare models. Can be one of the built-in methods or a
  function to compute the returned metrics (see Details).

- ...:

  models to be compared. Must be models created by serosv. If models'
  names are not provided, indices will be used instead for the \`model\`
  column in the returned data.frame.

## Value

a data.frame of 4 columns

- label:

  name or index of the model

- type:

  model type of the given model (a serosv model name)

- AIC:

  AIC value for the model (lower value indicates better fit)

- BIC:

  BIC value for the model (lower value indicates better fit)

## Details

Built-in comparison methods include: - computing AIC and BIC, which
returns AIC, BIC values of the model if available - cross validation,
which reutns
