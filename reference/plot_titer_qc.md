# Quality control plot

Visualize for each sample, whether titer estimates can be computed at
the tested dilutions.

## Usage

``` r
plot_titer_qc(x, n_plates = 18, n_samples = 22, n_dilutions = 3)
```

## Arguments

- x:

  \- output of \`to_titer()\`

- n_plates:

  \- maximum number of plates to plot

- n_samples:

  \- maximum number of samples per plate to plot

- n_dilutions:

  \- number of dilutions used for testing

## Details

Each sample is represented by a \`n_estimates x n_dilutions\` grid where
cell color indicate estimate availability (green = estimate available,
orange = result too low, red = result too high)

These sample grids are arranged in columns where each column represent
samples from a plate

The figure below demonstrates the interpretation of the plot.
![](figures/interpret_titer_qc.png)
