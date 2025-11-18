# A local polynomial model.

Refers to section 7.1. and 7.2.

## Usage

``` r
lp_model(data, kern = "tcub", nn = 0, h = 0, deg = 2)
```

## Arguments

- data:

  the input data frame, must either have \`age\`, \`pos\`, \`tot\`
  columns (for aggregated data) OR \`age\`, \`status\` for (linelisting
  data)

- kern:

  Weight function, default = "tcub". Other choices are "rect", "trwt",
  "tria", "epan", "bisq" and "gauss". Choices may be restricted when
  derivatives are required; e.g. for confidence bands and some bandwidth
  selectors.

- nn:

  Nearest neighbor component of the smoothing parameter. Default value
  is 0.7, unless either h is provided, in which case the default is 0.

- h:

  The constant component of the smoothing parameter. Default: 0.

- deg:

  Degree of polynomial to use. Default: 2.

## Value

a list of class lp_model with 6 items

- datatype:

  type of datatype used for model fitting (aggregated or linelisting)

- df:

  the dataframe used for fitting the model

- pi:

  fitted locfit object for pi

- eta:

  fitted locfit object for eta

- sp:

  seroprevalence

- foi:

  force of infection

## See also

\[locfit::locfit()\] for more information on the fitted locfit object

## Examples

``` r
df <- mumps_uk_1986_1987
model <- lp_model(
  df,
  nn=0.7, kern="tcub"
  )
plot(model)

```
