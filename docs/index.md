# serosv

`serosv` is an easy-to-use and efficient tool to estimate infectious
diseases parameters (seroprevalence and force of infection) using
serological data. The current version is mostly based on the book
“Modeling Infectious Disease Parameters Based on Serological and Social
Contact Data – A Modern Statistical Perspective” by [Hens et al., 2012
Springer](https://link.springer.com/book/10.1007/978-1-4614-4072-7).

## Installation

You can install the development version of serosv with:

``` r
# install.packages("devtools")
devtools::install_github("OUCRU-Modelling/serosv")
```

## Feature overview

### Datasets

`serosv` contains 15 built-in serological datasets as provided by [Hens
et al., 2012
Springer](https://link.springer.com/book/10.1007/978-1-4614-4072-7).
Simply call the name to load a dataset, for example:

``` r
rubella <- rubella_uk_1986_1987
```

### Methods

The following methods are available to estimate seroprevalence and force
of infection.

Parametric approaches:

- Frequentist methods:

  - Polynomial models:

    - Muench’s model
    - Griffiths’ model
    - Grenfell and Anderson’s model

  - Nonlinear models:

    - Farrington’s model
    - Weibull model

  - Fractional polynomial models

- Bayesian methods:

  - Hierarchical Farrington model

  - Hierarchical log-logistic model

Nonparametric approaches:

- Local estimation by polynomials

Semiparametric approaches:

- Penalized splines:

  - Penalized likelihood framework

  - Generalized Linear Mixed Model framework

## Demo

### Fitting rubella data from the UK

Load the rubella in UK dataset.

``` r
library(serosv)
rubella <- rubella_uk_1986_1987
```

Find the power for the best second degree fractional polynomial with
monotonicity constraint and a logit link function. The power appears to
be (-0.9,-0.9).

``` r
rubella_mod <- fp_model(
  rubella,
  p=list(
    p_range=seq(-2,3,0.1),
    degree=2
  ), 
  monotonic = T, link="logit"
)
rubella_mod
#> Fractional polynomial model 
#> 
#> Input type:  aggregated 
#> Powers:  -0.9, -0.9 
#> 
#> Call:  glm(formula = as.formula(formulate(curr_p)), family = binomial(link = link))
#> 
#> Coefficients:
#>               (Intercept)                I(age^-0.9)  
#>                     4.342                     -4.696  
#> I(I(age^-0.9) * log(age))  
#>                    -9.845  
#> 
#> Degrees of Freedom: 43 Total (i.e. Null);  41 Residual
#> Null Deviance:       1369 
#> Residual Deviance: 37.58     AIC: 210.1
```

Visualize the model

``` r
plot(rubella_mod)
```

![](reference/figures/README-unnamed-chunk-5-1.png)

### Fitting Parvo B19 data from Finland

``` r
library(dplyr)
#> Warning: package 'dplyr' was built under R version 4.3.1
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
parvob19 <- parvob19_fi_1997_1998

# for linelisting data, either transform it to aggregated
transform_data(
  parvob19$age, 
  parvob19$seropositive,
  stratum_col = "age") |>
  polynomial_model(k = 1) |>
  plot()
```

![](reference/figures/README-unnamed-chunk-6-1.png)

``` r

# or fit data as is
parvob19 |>
  rename(status = seropositive) |>
  polynomial_model(k = 1) |>
  plot()
```

![](reference/figures/README-unnamed-chunk-6-2.png)
