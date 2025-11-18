# Plot output for age_time_model

Plot output for age_time_model

## Usage

``` r
# S3 method for class 'age_time_model'
plot(x, ...)
```

## Arguments

- x:

  \- a \`age_time_model\` object

- ...:

  arbitrary params. Supported options include:

  - `facet`: Whether to facet the plot by group.

  - `modtype`: Which model to plot, either `"monotonized"` or
    `"non-monotonized"`.

  - `le`: Number of bins used to generate the x-axis; higher values
    produce smoother curves.

  - `cex`: Adjusts the size of data points (only when `facet = TRUE`).

## Value

ggplot object
