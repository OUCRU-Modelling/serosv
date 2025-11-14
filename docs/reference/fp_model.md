# A fractional polynomial model.

Refers to section 6.2.

## Usage

``` r
fp_model(data, p, link = "logit")
```

## Arguments

- data:

  the input data frame, must either have \`age\`, \`pos\`, \`tot\`
  columns (for aggregated data) OR \`age\`, \`status\` for (linelisting
  data)

- p:

  the powers of the predictor.

- link:

  the link function for model. Defaulted to "logit".

## Value

a list of class fp_model with 5 items

- datatype:

  type of data used for fitting model (aggregated or linelisting)

- df:

  the dataframe used for fitting the model

- info:

  a fitted glm model

- sp:

  seroprevalence

- foi:

  force of infection

## See also

\[stats::glm()\] for more information on glm object

## Examples

``` r
df <- hav_be_1993_1994
model <- fp_model(
  df,
  p=c(1.5, 1.6), link="cloglog")
plot(model)

```
