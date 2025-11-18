# The Weibull model.

Refers to section 6.1.2.

## Usage

``` r
weibull_model(data)
```

## Arguments

- data:

  the input data frame, must either have \`t\`, \`pos\`, \`tot\` column
  for aggregated data OR \`t\`, \`status\` for linelisting data

## Value

list of class weibull_model with the following items

- datatype:

  type of datatype used for model fitting (aggregated or linelisting)

- df:

  the dataframe used for fitting the model

- info:

  fitted "glm" object

- sp:

  seroprevalence

- foi:

  force of infection

## See also

\[stats::glm()\] for more information on the fitted "glm" object

## Examples

``` r
df <- hcv_be_2006[order(hcv_be_2006$dur), ]
df$t <- df$dur
df$status <- df$seropositive
model <- weibull_model(df)
plot(model)

```
