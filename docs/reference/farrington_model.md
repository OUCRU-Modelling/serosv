# The Farrington (1990) model.

Refers to section 6.1.2.

## Usage

``` r
farrington_model(data, start, fixed = list())
```

## Arguments

- data:

  the input data frame, must either have \`age\`, \`pos\`, \`tot\`
  columns (for aggregated data) OR \`age\`, \`status\` for (linelisting
  data)

- start:

  Named list of vectors or single vector. Initial values for optimizer.

- fixed:

  Named list of vectors or single vector. Parameter values to keep fixed
  during optimization.

## Value

a list of class farrington_model with 5 items

- datatype:

  type of datatype used for model fitting (aggregated or linelisting)

- df:

  the dataframe used for fitting the model

- info:

  fitted "mle" object

- sp:

  seroprevalence

- foi:

  force of infection

## See also

\[stats4::mle()\] for more information on the fitted mle object

## Examples

``` r
df <- rubella_uk_1986_1987
model <- farrington_model(
  df,
  start=list(alpha=0.07,beta=0.1,gamma=0.03)
  )
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
plot(model)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.

```
