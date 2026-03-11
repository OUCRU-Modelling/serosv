# Model selection

``` r
library(serosv)
library(dplyr)
library(magrittr)
```

``` r
data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age), ] 

aggregated <- transform_data(data$age, data$seropositive, stratum_col = "age")
```

## Generate models comparison `data.frame`

Function
[`compare_models()`](https://oucru-modelling.github.io/serosv/reference/compare_models.md)
is used for quickly computing comparison metrics for a set of models on
a given dataset.

The function takes the following arguments:

- `data` the dataset that will be fitted to the models

- `method` the method to generate comparison metrics. It can be the name
  of one of the built-in methods, or a user-defined function.

- `…` functions that take a data and return a fitted `serosv` model.

It will then return a `data.frame` with the following columns

- `label` model identifier. Either user defined name or index based on
  the order provided.

- `type` type of model (a `serosv` model class)

- metrics depending on the method selected

**Built-in `method`**

`serosv` currently provide 2 built-in metrics generating methods

- `"AIC/BIC"` which returns fitted model’s AIC, BIC and Log-likelihood
  (where applicable)

- `"CV"` which split the data into “train” and “test” set then return
  logloss and AUC from model’s prediction on the test set

**Sample usage**

``` r
# ----- Return AIC, BIC ----- 
aic_bic_out <- compare_models(
  data = data %>% rename(status=seropositive),
  method = "AIC/BIC",
  griffith = ~polynomial_model(.x, k=3),
  penalized_spline = penalized_spline_model,
  farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
  local_polynomial = lp_model # expect to not return any values
  ) %>% suppressWarnings()

# ----- Return Cross-validation metrics ----- 
cv_out <- compare_models(
    data %>% rename(status=seropositive),
    method = "CV",
    griffith = ~polynomial_model(.x, k=3),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
    local_polynomial = lp_model 
  ) %>% suppressWarnings()

aic_bic_out
```

    ## # A tibble: 4 × 7
    ##   label            type                     AIC   BIC logLik    df mod_out   
    ##   <chr>            <chr>                  <dbl> <dbl>  <dbl> <dbl> <list>    
    ## 1 griffith         polynomial_model       1347. 1362.  -671.  3    <plynml_m>
    ## 2 penalized_spline penalized_spline_model 1330. 1365.  -658.  6.88 <pnlzd_s_>
    ## 3 farrington       farrington_model       1337.   NA   -665.  3    <frrngtn_>
    ## 4 local_polynomial lp_model                 NA    NA     NA  NA    <lp_model>

``` r
cv_out
```

    ## # A tibble: 4 × 5
    ##   label            type                   logloss   auc mod_out   
    ##   <chr>            <chr>                    <dbl> <dbl> <list>    
    ## 1 griffith         polynomial_model         -209. 0.700 <plynml_m>
    ## 2 penalized_spline penalized_spline_model   -132. 0.673 <pnlzd_s_>
    ## 3 farrington       farrington_model         -133. 0.700 <frrngtn_>
    ## 4 local_polynomial lp_model                 -131. 0.700 <lp_model>

With aggregated data

``` r
# ----- Return AIC, BIC ----- 
aic_bic_out <- compare_models(
  data = aggregated,
  method = "AIC/BIC",
  griffith = ~polynomial_model(.x, k=3),
  penalized_spline = penalized_spline_model,
  farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
  local_polynomial = lp_model # expect to not return any values
  ) %>% suppressWarnings()

# ----- Return Cross-validation metrics ----- 
cv_out <- compare_models(
    data = aggregated,
    method = "CV",
    griffith = ~polynomial_model(.x, k=3),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
    local_polynomial = lp_model 
  ) %>% suppressWarnings()


aic_bic_out
```

    ## # A tibble: 4 × 7
    ##   label            type                     AIC   BIC logLik    df mod_out   
    ##   <chr>            <chr>                  <dbl> <dbl>  <dbl> <dbl> <list>    
    ## 1 griffith         polynomial_model        484.  494.  -239.  3    <plynml_m>
    ## 2 penalized_spline penalized_spline_model  257.  271.  -124.  4.37 <pnlzd_s_>
    ## 3 farrington       farrington_model       1337.   NA   -665.  3    <frrngtn_>
    ## 4 local_polynomial lp_model                 NA    NA     NA  NA    <lp_model>

``` r
cv_out
```

    ## # A tibble: 4 × 5
    ##   label            type                      mae logloss mod_out   
    ##   <chr>            <chr>                   <dbl>   <dbl> <list>    
    ## 1 griffith         polynomial_model       1.82     -85.0 <plynml_m>
    ## 2 penalized_spline penalized_spline_model 0.0261   -43.2 <pnlzd_s_>
    ## 3 farrington       farrington_model       0.234    -42.1 <frrngtn_>
    ## 4 local_polynomial lp_model               0.0161   -42.7 <lp_model>
