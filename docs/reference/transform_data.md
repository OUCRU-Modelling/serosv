# Aggregate data

Generate a dataframe with \`t\`, \`pos\` and \`tot\` columns from \`t\`
and \`seropositive\` vectors.

## Usage

``` r
transform_data(data, stratum_col = "age", status_col = "status")
```

## Arguments

- data:

  a data frame with columns for age and serostatus

- stratum_col:

  name of the column to stratify by (default to "age")

- status_col:

  name of the column for serostatus

## Value

dataframe in aggregated format

## Examples

``` r
df <- hcv_be_2006
hcv_df <- transform_data(df, stratum_col="dur", status_col="seropositive")
hcv_df
#> # A tibble: 116 × 3
#>      age   pos   tot
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
