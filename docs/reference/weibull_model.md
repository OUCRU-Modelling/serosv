# The Weibull model.

Model seroprevalence as a function of age using the Weibull model, where
the force of infection is assumed to vary monotonically with age.

## Usage

``` r
weibull_model(
  data,
  age_col = "age",
  pos_col = "pos",
  tot_col = "tot",
  status_col = "status",
  ...
)
```

## Arguments

- data:

  the input data frame, must either have columns for \`age\`, \`pos\`,
  \`tot\` (for aggregated data) OR \`age\`, \`status\` (for linelisting
  data)

- age_col:

  name of the \`age\` column (default age_col="age")

- pos_col:

  name of the \`pos\` column (default pos_col="pos")

- tot_col:

  name of the \`tot\` column (default tot_col="tot")

- status_col:

  name of the \`status\` column (default status_col="status")

- ...:

  additional arguments to be passed to \`glm()\` function that fits the
  model

## Value

list of class weibull_model with the following items

- datatype:

  type of datatype used for model fitting (aggregated or linelisting)

- df:

  the dataframe used for fitting the model

- info:

  fitted "glm" object

- sp:

  estimated seroprevalence

- foi:

  estimated force of infection

- foi_mod:

  function to generate FoI given age, and parameter values

## Details

For a Weibull model, the prevalence is given by \$\$ \pi (a) = 1 - e^{ -
\beta_0 a ^ {\beta_1}} \$\$ Where \\a\\ is the age, which may refer to
biological age or a time scale of interest (e.g., time since
vaccination).

Which implies the force of infection to be the monotonic function \$\$
\lambda(a) = \beta_0 \beta_1 a^{\beta_1 - 1} \$\$

Refer to section 6.1.2. of the the book by Hens et al. (2012) for
further details.

## References

Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme, and
Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
Serological and Social Contact Data: A Modern Statistical Perspective.
tatistics for Biology and Health. Springer New York.
[doi:10.1007/978-1-4614-4072-7](https://doi.org/10.1007/978-1-4614-4072-7)
.

## See also

\[stats::glm()\] for more information on the fitted "glm" object

## Examples

``` r
df <- hcv_be_2006[order(hcv_be_2006$dur), ]
model <- weibull_model(df, age_col="dur", status_col="seropositive")
plot(model)

```
