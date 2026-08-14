# Generate table of metrics for model comparison

Generate table of metrics for model comparison

## Usage

``` r
compare_models(data, method = "AIC/BIC", method_args = list(), ...)
```

## Arguments

- data:

  input data to fit into the models

- method:

  method to compare models. Can be one of the built-in methods or a
  function to compute the returned metrics (see Details).

- method_args:

  additional arguments to be passed to the method function.

- ...:

  models to be compared. Must be models created by serosv. If models'
  names are not provided, indices will be used instead for the \`model\`
  column in the returned data.frame.

## Value

a data.frame with the following columns

- label:

  name or index of the model

- type:

  model type of the given model (a serosv model name)

- mod_out:

  the fitted models

- plots:

  the plots for each of the fitted model

- metrics columns:

  the columns for metrics of comparison, the number of which depends on
  the function that generate these metrics

## Details

Built-in comparison methods include:

- computing AIC and BIC, which returns AIC, BIC values of the model if
  available

- cross validation (perform k-fold validation), which returns MSE and
  logloss (negative log Binomial likelihood) for aggregated data, or AUC
  and logloss (negative log Bernoulli likelihood) for linelisting data

## Examples

``` r
comparison_table <- suppressWarnings(
  compare_models(
    data = hav_bg_1964,
    method = "CV",
    polynomial_mod = ~polynomial_model(.x, k=1),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, start=list(alpha=0.3,beta=0.1,gamma=0.03))
  )
)
#> Error in map2(.x, vec_index(.x), .f, ...): ℹ In index: 3.
#> ℹ With name: farrington.
#> Caused by error in `mutate()`:
#> ℹ In argument: `plots = list(plot(out) + ggtitle(paste("Fitted model
#>   using", class(out))))`.
#> Caused by error in `quantile.default()`:
#> ! missing values and NaN's not allowed if 'na.rm' is FALSE
# view table of metrics
comparison_table
#> Error: object 'comparison_table' not found
# view the model fitted with the whole dataset
comparison_table$plots
#> Error: object 'comparison_table' not found
```
