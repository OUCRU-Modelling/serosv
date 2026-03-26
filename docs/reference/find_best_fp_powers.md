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

## Value

list of 3 elements:

- p:

  The best power for fp model.

- deviance:

  Deviance of the best fitted model.

- model:

  The best model fitted
