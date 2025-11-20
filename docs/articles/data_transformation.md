# Data transformation

``` r
library(serosv)
```

## Aggregate data

Use function
[`transform_data()`](https://oucru-modelling.github.io/serosv/reference/transform_data.md)
to convert line listing data (data with `t`, `seropositive` columns) to
aggregated data (data with `t`, `pos`, `tot` columns)

Arguments:

- `t` the time vector

- `spos` the seropositive vector (TRUE/FALSE or 1/0)

- `stratum_col` new name for the time vector, default to
  `stratum_col = "t"`

``` r
linelisting_data <- hcv_be_2006
aggregated_data <- transform_data(linelisting_data$dur, linelisting_data$seropositive)
aggregated_data
#> # A tibble: 116 × 3
#>        t   pos   tot
#>    <dbl> <int> <int>
#>  1   0.1     0     1
#>  2   0.3     0     1
#>  3   0.4     0     3
#>  4   0.5     0     1
#>  5   0.7     2     5
#>  6   0.8     3     4
#>  7   0.9     3     4
#>  8   1       5     7
#>  9   1.1     1     2
#> 10   1.2     3     4
#> # ℹ 106 more rows
```
