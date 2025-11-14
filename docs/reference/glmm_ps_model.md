# Fitting penalized spline model with Generalized Linear Mixture Model framework

Fitting penalized spline model with Generalized Linear Mixture Model
framework

## Usage

``` r
glmm_ps_model(age, res, s = "bs", link = "logit")
```

## Arguments

- age:

  the age vector.

- res:

  the serostatus vector.

- s:

  smoothing basis to use.

- link:

  link function to use.

## Examples

``` r
data <- parvob19_be_2001_2003
mod <- glmm_ps_model(data$age, data$seropositive)
#> 
#>  Maximum number of PQL iterations:  20 
#> iteration 1
#> iteration 2
#> iteration 3
#> iteration 4
mod$gam$info
#> NULL
```
