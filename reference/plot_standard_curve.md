# Visualize standard curves for each plate

Visualize standard curves for each plate

## Usage

``` r
plot_standard_curve(
  x,
  facet = TRUE,
  xlab = "log10(concentration)",
  ylab = "Optical density",
  datapoint_size = 2
)
```

## Arguments

- x:

  \- output of \`to_titer()\`

- facet:

  \- whether to faceted by plates or plot all standard curves on a
  single plot

- xlab:

  \- label of the x axis

- ylab:

  \- label of the y axis

- datapoint_size:

  \- size of the data point (only applicable when \`facet=TRUE\`)
