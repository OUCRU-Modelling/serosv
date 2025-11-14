# Compute confidence interval for time age model

Compute confidence interval for time age model

## Usage

``` r
compute_ci.age_time_model(x, ci = 0.95, le = 100, ...)
```

## Arguments

- x:

  \- serosv models

- ci:

  \- confidence interval

- le:

  \- number of data for computing confidence interval

- ...:

  \- arbitrary argument

## Value

confidence interval dataframe with n_group x 3 cols, the columns are
\`group\`, \`sp_df\`, \`foi_df\`
