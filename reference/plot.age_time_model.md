# Plot output for age_time_model

Plot output for age_time_model

## Usage

``` r
# S3 method for class 'age_time_model'
plot(x, cex = 10, le = 100, facet = TRUE, ...)
```

## Arguments

- x:

  a \`age_time_model\` object

- cex:

  adjust size of the datapoints (only when `facet = TRUE`)

- le:

  number of bins used to generate the x-axis; higher values produce
  smoother curves

- facet:

  whether to facet the plot by group

- ...:

  arbitrary params

- modtype:

  specify which model to plot, either `"monotonized"` or
  `"non-monotonized"`

## Value

ggplot object
