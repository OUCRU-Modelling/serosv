# Returns the powers of the fractional polynomial model which has the lowest deviance score.

Return the best powers for a given degree

## Usage

``` r
find_best_fp_powers(data, p, mc, degree, link = "logit")
```

## Arguments

- data:

  the input data frame, must either have columns for \`age\`, \`pos\`,
  \`tot\` (for aggregated data) OR \`age\`, \`status\` (for linelisting
  data)

- p:

  a powers sequence to be tested.

- mc:

  indicates if the returned model should be monotonic.

- degree:

  the maximum degree (i.e. number of power terms) to search for the best
  model. Recommended to be \<= 2.

- link:

  the link function. Defaulted to "logit".

- age_col:

  name of the \`age\` column (default age_col="age").

- pos_col:

  name of the \`pos\` column (default pos_col="pos").

- tot_col:

  name of the \`tot\` column (default tot_col="tot").

- status_col:

  name of the \`status\` column (default status_col="status").

## Value

list of 3 elements:

- p:

  The best power for fp model.

- deviance:

  Deviance of the best fitted model.

- model:

  The best model fitted

## Examples

``` r
df <- hav_be_1993_1994
best_p <- find_best_fp_powers(
df,
p=seq(-2,3,0.1), mc=FALSE, degree=2, link="cloglog"
)
#> Error in find_best_fp_powers(df, p = seq(-2, 3, 0.1), mc = FALSE, degree = 2,     link = "cloglog"): could not find function "find_best_fp_powers"
best_p
#> Error: object 'best_p' not found
```
