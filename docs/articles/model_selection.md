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
  griffith = ~polynomial_model(.x, k=2),
  penalized_spline = penalized_spline_model,
  farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
  local_polynomial = lp_model # expect to not return any values
  ) %>% suppressWarnings()

# ----- Return Cross-validation metrics ----- 
cv_out <- compare_models(
    data %>% rename(status=seropositive),
    method = "CV",
    griffith = ~polynomial_model(.x, k=2),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
    local_polynomial = lp_model 
  ) %>% suppressWarnings()

aic_bic_out
```

    ## # A tibble: 4 × 7
    ##   label            type                     AIC   BIC logLik    df mod_out   
    ##   <chr>            <chr>                  <dbl> <dbl>  <dbl> <dbl> <list>    
    ## 1 griffith         polynomial_model       1353. 1363.  -674.  2    <plynml_m>
    ## 2 penalized_spline penalized_spline_model 1330. 1365.  -658.  6.88 <pnlzd_s_>
    ## 3 farrington       farrington_model       1337.   NA   -665.  3    <frrngtn_>
    ## 4 local_polynomial lp_model                 NA    NA     NA  NA    <lp_model>

``` r
cv_out
```

    ## # A tibble: 4 × 4
    ##   label            logloss   auc type                  
    ##   <chr>              <dbl> <dbl> <chr>                 
    ## 1 griffith            277. 0.594 polynomial_model      
    ## 2 penalized_spline    197. 0.601 penalized_spline_model
    ## 3 farrington          356. 0.600 farrington_model      
    ## 4 local_polynomial    237. 0.588 lp_model

With aggregated data

``` r
# ----- Return AIC, BIC ----- 
aic_bic_out <- compare_models(
  data = aggregated,
  method = "AIC/BIC",
  griffith = ~polynomial_model(.x, k=2),
  penalized_spline = penalized_spline_model,
  farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
  local_polynomial = lp_model # expect to not return any values
  ) %>% suppressWarnings()

# ----- Return Cross-validation metrics (default with 4 folds) ----- 
cv_out <- compare_models(
    data = aggregated,
    method = "CV",
    griffith = ~polynomial_model(.x, k=2),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
    local_polynomial = lp_model 
  ) %>% suppressWarnings()


aic_bic_out
```

    ## # A tibble: 4 × 7
    ##   label            type                     AIC   BIC logLik    df mod_out   
    ##   <chr>            <chr>                  <dbl> <dbl>  <dbl> <dbl> <list>    
    ## 1 griffith         polynomial_model        489.  496.  -243.  2    <plynml_m>
    ## 2 penalized_spline penalized_spline_model  257.  271.  -124.  4.37 <pnlzd_s_>
    ## 3 farrington       farrington_model       1337.   NA   -665.  3    <frrngtn_>
    ## 4 local_polynomial lp_model                 NA    NA     NA  NA    <lp_model>

``` r
cv_out
```

    ## # A tibble: 4 × 4
    ##   label              mse logloss type                  
    ##   <chr>            <dbl>   <dbl> <chr>                 
    ## 1 griffith         0.217   178.  polynomial_model      
    ## 2 penalized_spline 0.123    79.7 penalized_spline_model
    ## 3 farrington       0.117    65.5 farrington_model      
    ## 4 local_polynomial 0.117    71.3 lp_model
