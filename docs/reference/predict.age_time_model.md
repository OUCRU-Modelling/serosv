# Predict from the age_time_mdoel

Predict from the age_time_mdoel

## Usage

``` r
# S3 method for class 'age_time_model'
predict(object, newdata, modtype = "monotonized", ...)
```

## Arguments

- object:

  serosv models

- newdata:

  data.frame with age column to generate prediction

- modtype:

  either "monotonized" (to predict using monotonized model) or
  "non-monotonized"

- ...:

  arbitrary argument

## Value

confidence interval dataframe with n_group x 3 cols, the columns are
\`group\`, \`sp_df\`, \`foi_df\`
