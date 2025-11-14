# Penalized Spline model

Penalized Spline model

## Usage

``` r
penalized_spline_model(
  data,
  s = "bs",
  link = "logit",
  framework = "pl",
  sp = NULL
)
```

## Arguments

- data:

  the input data frame, must either have \`age\`, \`pos\`, \`tot\`
  column for aggregated data OR \`age\`, \`status\` for linelisting data

- s:

  smoothing basis to use

- link:

  link function to use

- framework:

  which approach to fit the model ("pl" for penalized likelihood
  framework, "glmm" for generalized linear mixed model framework)

- sp:

  smoothing parameter

## Value

a list of class penalized_spline_model with 6 attributes

- datatype:

  type of datatype used for model fitting (aggregated or linelisting)

- df:

  the dataframe used for fitting the model

- framework:

  either pl or glmm

- info:

  fitted "gam" model when framework is pl or "gamm" model when framework
  is glmm

- sp:

  seroprevalence

- foi:

  force of infection

## See also

\[mgcv::gam()\], \[mgcv::gamm()\] for more information the fitted gam
and gamm model

## Examples

``` r
data <- parvob19_be_2001_2003
data$status <- data$seropositive
model <- penalized_spline_model(data, framework="glmm")
#> 
#>  Maximum number of PQL iterations:  20 
#> iteration 1
#> iteration 2
#> iteration 3
#> iteration 4
model$info$gam
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> spos ~ s(age, bs = s, sp = sp)
#> 
#> Estimated degrees of freedom:
#> 5.8  total = 6.8 
#> 
plot(model)
```
