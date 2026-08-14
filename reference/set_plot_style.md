# Helper to adjust styling of a plot

Helper to adjust styling of a plot

## Usage

``` r
set_plot_style(
  sero = "blueviolet",
  sero_ci = "royalblue1",
  foi = "#fc0328",
  foi_ci = "#fc0328",
  sero_line = "solid",
  foi_line = "dashed",
  xlabel = "Age"
)
```

## Arguments

- sero:

  color for seroprevalence line

- sero_ci:

  color for confidence intervals of seroprevalence

- foi:

  color for force of infection line

- foi_ci:

  color for confidence intervals of FoI

- sero_line:

  linetype for seroprevalence line

- foi_line:

  linetype for force of infection line

- xlabel:

  x label

## Value

list of updated aesthetic values
