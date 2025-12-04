# A fractional polynomial model.

Fractional polynomial model is a generalization of polynomial models
where the power of the terms can be fractions, allowing more flexibility
and better fit for data where asymptotic behavior is expected.

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

## Details

Instead of a polynomial, the linear predictor is now defined as \$\$
\eta_m(a, \beta, p_1, p_2, ...,p_m) = \Sigma^m\_{i=0} \beta_i H_i(a)
\$\$ Where \\(m\\) is an integer, \\(p_1 \le p_2 \le... \le p_m\\) is a
sequence of powers, and \\(H_i(a)\\) is a transformation given by

\$\$ H_i = \begin{cases} a^{p_i} & \text{ if } p_i \neq p\_{i-1}, \\\\
H\_{i-1}(a) \times log(a) & \text{ if } p_i = p\_{i-1}, \end{cases} \$\$

Refers to section 6.2. of the the book by Hens et al. (2012) for further
details.

## References

Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme, and
Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
Serological and Social Contact Data: A Modern Statistical Perspective.
tatistics for Biology and Health. Springer New York.
[doi:10.1007/978-1-4614-4072-7](https://doi.org/10.1007/978-1-4614-4072-7)
.

## See also

\[stats::glm()\] for more information on glm object

\[polynomial_models()\]

## Examples

``` r
df <- hav_be_1993_1994
model <- fp_model(
  df,
  p=c(1.5, 1.6), link="cloglog")
plot(model)

```
