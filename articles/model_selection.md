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

- `method_args` the list of additional arguments for method function.

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
  logloss and AUC from model’s prediction on the test set. Number of
  folds can be configured via `k` argument (see **Sample usage** below).

**Sample usage**

``` r
# ----- Return AIC, BIC ----- 
aic_bic_out <- compare_models(
  data = data,
  method = "AIC/BIC",
  griffith = ~polynomial_model(.x, k=2),
  penalized_spline = penalized_spline_model,
  farrington = ~farrington_model(.x, 
                                 start=list(beta=0.12,gamma=0.05),
                                 fixed=list(alpha=0.1)),
  local_polynomial = lp_model # expect to not return any values
  ) %>% suppressWarnings()

# ----- Return Cross-validation metrics ----- 
cv_out <- compare_models(
    data,
    method = "CV",
    griffith = ~polynomial_model(.x, k=2),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, 
                                   start=list(beta=0.12,gamma=0.05),
                                   fixed=list(alpha=0.1)),
    local_polynomial = lp_model,
    method_args = list(
      k = 4 # adjust the fold k for the CV function
    )
  ) %>% suppressWarnings()

aic_bic_out
```

    ## # A tibble: 4 × 8
    ##   label            type             AIC   BIC logLik    df mod_out    plots     
    ##   <chr>            <chr>          <dbl> <dbl>  <dbl> <dbl> <list>     <list>    
    ## 1 griffith         polynomial_mo… 1353. 1363.  -674.  2    <plynml_m> <ggplt2::>
    ## 2 penalized_spline penalized_spl… 1330. 1365.  -658.  6.88 <pnlzd_s_> <ggplt2::>
    ## 3 farrington       farrington_mo… 1358.   NA   -677.  2    <frrngtn_> <ggplt2::>
    ## 4 local_polynomial lp_model         NA    NA     NA  NA    <lp_model> <ggplt2::>

``` r
cv_out
```

    ## # A tibble: 4 × 6
    ##   label            logloss   auc type                   mod_out    plots     
    ##   <chr>              <dbl> <dbl> <chr>                  <list>     <list>    
    ## 1 griffith            169. 0.701 polynomial_model       <plynml_m> <ggplt2::>
    ## 2 penalized_spline    167. 0.684 penalized_spline_model <pnlzd_s_> <ggplt2::>
    ## 3 farrington          170. 0.702 farrington_model       <frrngtn_> <ggplt2::>
    ## 4 local_polynomial    165. 0.702 lp_model               <lp_model> <ggplt2::>

With aggregated data

``` r
# ----- Return AIC, BIC ----- 
aic_bic_out <- compare_models(
  data = aggregated,
  method = "AIC/BIC",
  griffith = ~polynomial_model(.x, k=2),
  penalized_spline = penalized_spline_model,
  farrington = ~farrington_model(.x, 
                                 start=list(beta=0.1,gamma=0.03),
                                 fixed=list(alpha=0.07)),
  local_polynomial = lp_model # expect to not return any values
  ) %>% suppressWarnings()

# ----- Return Cross-validation metrics (default with 4 folds) ----- 
cv_out <- compare_models(
    data = aggregated,
    method = "CV",
    griffith = ~polynomial_model(.x, k=2),
    penalized_spline = penalized_spline_model,
    farrington = ~farrington_model(.x, 
                                 start=list(beta=0.1,gamma=0.03),
                                 fixed=list(alpha=0.07)),
    local_polynomial = lp_model,
    method_args = list(
      k = 4 # adjust the fold k for the CV function
    )
  ) %>% suppressWarnings()

aic_bic_out
```

    ## # A tibble: 4 × 8
    ##   label            type             AIC   BIC logLik    df mod_out    plots     
    ##   <chr>            <chr>          <dbl> <dbl>  <dbl> <dbl> <list>     <list>    
    ## 1 griffith         polynomial_mo…  489.  496.  -243.  2    <plynml_m> <ggplt2::>
    ## 2 penalized_spline penalized_spl…  467.  490.  -227.  6.88 <pnlzd_s_> <ggplt2::>
    ## 3 farrington       farrington_mo…  486.   NA   -241.  2    <frrngtn_> <ggplt2::>
    ## 4 local_polynomial lp_model         NA    NA     NA  NA    <lp_model> <ggplt2::>

``` r
cv_out
```

    ## # A tibble: 4 × 6
    ##   label              mse logloss type                   mod_out    plots     
    ##   <chr>            <dbl>   <dbl> <chr>                  <list>     <list>    
    ## 1 griffith         0.110    61.7 polynomial_model       <plynml_m> <ggplt2::>
    ## 2 penalized_spline 0.110    58.4 penalized_spline_model <pnlzd_s_> <ggplt2::>
    ## 3 farrington       0.107    60.5 farrington_model       <frrngtn_> <ggplt2::>
    ## 4 local_polynomial 0.109    60.2 lp_model               <lp_model> <ggplt2::>

## Generate custom metrics

The users can also provide a custom function to generate the comparison
metrics.

This function must accepts 2 parameters:

- `dat` the input data

- `mod_func` a function that takes an input data and returns a `serosv`
  model

And it must returns a data frame with 1 row where each column represents
one metric.

***Example:***

The following implements holdout validation and returns MAE:

``` r
generate_mae <- function(dat, mod_func){
  n_train <- round(nrow(dat)*0.8)
  train <- dat[1:n_train,]
  test <- dat[n_train:nrow(dat),]
  
  fit <- mod_func(dat)
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
    ## 1         griffith 0.3774195
    ## 2 penalized_spline 0.3249623
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
    ## 1         griffith 0.4603812
    ## 2 penalized_spline 0.4540285
    ## 3       farrington 0.4636040
    ## 4 local_polynomial 0.4536016
