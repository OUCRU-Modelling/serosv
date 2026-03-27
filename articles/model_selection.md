# Model selection

``` r
library(serosv)
library(dplyr)
library(magrittr)
```

``` r
data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age), ] %>% 
  rename(status = seropositive) 

aggregated <- transform_data(data, stratum_col = "age", status_col="status")
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
  data = data,
  method = "AIC/BIC",
  griffith = ~polynomial_model(.x, k=2),
  penalized_spline = penalized_spline_model,
  farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
  local_polynomial = lp_model # expect to not return any values
  ) %>% suppressWarnings()

# ----- Return Cross-validation metrics ----- 
cv_out <- compare_models(
    data,
    method = "CV",
    griffith = ~polynomial_model(.x, k=2),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
    local_polynomial = lp_model 
  ) %>% suppressWarnings()

aic_bic_out
```

    ## # A tibble: 4 × 8
    ##   label            type                  AIC   BIC logLik    df mod_out    plots
    ##   <chr>            <chr>               <dbl> <dbl>  <dbl> <dbl> <list>     <lis>
    ## 1 griffith         polynomial_model    1353. 1363.  -674.  2    <plynml_m> <gg> 
    ## 2 penalized_spline penalized_spline_m… 1330. 1365.  -658.  6.88 <pnlzd_s_> <gg> 
    ## 3 farrington       farrington_model    1337.   NA   -665.  3    <frrngtn_> <gg> 
    ## 4 local_polynomial lp_model              NA    NA     NA  NA    <lp_model> <gg>

``` r
cv_out
```

    ## # A tibble: 4 × 6
    ##   label            logloss   auc type                   mod_out    plots 
    ##   <chr>              <dbl> <dbl> <chr>                  <list>     <list>
    ## 1 griffith            277. 0.594 polynomial_model       <plynml_m> <gg>  
    ## 2 penalized_spline    197. 0.601 penalized_spline_model <pnlzd_s_> <gg>  
    ## 3 farrington          356. 0.600 farrington_model       <frrngtn_> <gg>  
    ## 4 local_polynomial    237. 0.588 lp_model               <lp_model> <gg>

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

    ## # A tibble: 4 × 8
    ##   label            type                  AIC   BIC logLik    df mod_out    plots
    ##   <chr>            <chr>               <dbl> <dbl>  <dbl> <dbl> <list>     <lis>
    ## 1 griffith         polynomial_model     489.  496.  -243.  2    <plynml_m> <gg> 
    ## 2 penalized_spline penalized_spline_m…  257.  271.  -124.  4.37 <pnlzd_s_> <gg> 
    ## 3 farrington       farrington_model    1337.   NA   -665.  3    <frrngtn_> <gg> 
    ## 4 local_polynomial lp_model              NA    NA     NA  NA    <lp_model> <gg>

``` r
cv_out
```

    ## # A tibble: 4 × 6
    ##   label              mse logloss type                   mod_out    plots 
    ##   <chr>            <dbl>   <dbl> <chr>                  <list>     <list>
    ## 1 griffith         0.217   178.  polynomial_model       <plynml_m> <gg>  
    ## 2 penalized_spline 0.123    79.7 penalized_spline_model <pnlzd_s_> <gg>  
    ## 3 farrington       0.117    65.5 farrington_model       <frrngtn_> <gg>  
    ## 4 local_polynomial 0.117    71.3 lp_model               <lp_model> <gg>

## Generate custom metrics

The users can also provide a custom function to generate the comparison
metrics.

This function must accepts 2 parameters:

- `dat` the input data

- `model_func` a function that takes an input data and returns a
  `serosv` model

And it must returns a data frame with 1 row where each column represents
one metric.

***Example:***

The following implements holdout validation and returns MAE:

``` r
generate_mae <- function(dat, model_func){
  n_train <- round(nrow(dat)*0.8)
  train <- dat[1:n_train,]
  test <- dat[n_train:nrow(dat),]
  
  fit <- model_func(dat)
  pred <- predict(fit, test[, 1, drop=FALSE])
  
  # handle error differently depending on datatype
  mae <- if(fit$datatype == "linelisting"){
    sum(abs(test$status - pred), na.rm=TRUE)/nrow(test)
  }else{
    sum(abs(test$pos/test$tot - pred), na.rm=TRUE)/nrow(test)
  }
  
  data.frame(
    mae = mae
  )
}
```

We can then run `compare_models` with the custom metrics function

``` r
compare_models(
    data = aggregated,
    method = generate_mae,
    griffith = ~polynomial_model(.x, k=2),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
    local_polynomial = lp_model 
  ) %>% suppressWarnings()
```

    ##              label       mae
    ## 1         griffith 1.9127037
    ## 2 penalized_spline 0.4801265
    ## 3       farrington 0.3710685
    ## 4 local_polynomial 0.3285755

``` r
compare_models(
    data = data,
    method = generate_mae,
    griffith = ~polynomial_model(.x, k=2),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, start=list(alpha=0.07,beta=0.1,gamma=0.03)),
    local_polynomial = lp_model 
  ) %>% suppressWarnings()
```

    ##              label       mae
    ## 1         griffith 1.6883332
    ## 2 penalized_spline 0.5215807
    ## 3       farrington 0.4636040
    ## 4 local_polynomial 0.4536016
