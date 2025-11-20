# Model selection

``` r
library(serosv)
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.3.1

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(magrittr)
```

``` r
data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age), ] %>% 
  rename(status = seropositive)
  
aggregated <- transform_data(data$age, data$status, stratum_col = "age")

# fit with linelisting data
model1 <- polynomial_model(data, type = "Muench")
# fit with aggregated data
model2 <- polynomial_model(aggregated, type = "Muench")
# fit with aggregated data
model3 <- polynomial_model(aggregated, type = "Griffith")
# fit with gam
model4 <- penalized_spline_model(aggregated)
```

## Generate models comparison `data.frame`

Function
[`compare_models()`](https://oucru-modelling.github.io/serosv/reference/compare_models.md)
is used for quickly computing AIC and BIC values for given model(s).

The function can take an arbitrary number of models and all models must
be created from `serosv`’s set of `*_model()` functions. It will then
return a `data.frame` of 4 columns:

- `model` model identifier. Either user defined name or index based on
  the order provided.

- `type` type of model (a `serosv` model class)

- `AIC` AIC value for the fitted model (if applicable)

- `BIC` AIC value for the fitted model (if applicable)

**Sample usage**

Compare 4 models defined above

``` r
# provide models with name
compare_models(muench_linelist = model1, muench_aggregated = model2, griffith = model3, penalized_spline = model4)
```

    ##               model                   type       AIC       BIC
    ## 1   muench_linelist       polynomial_model 1368.9096 1373.9280
    ## 2 muench_aggregated       polynomial_model  505.5269  508.8787
    ## 3          griffith       polynomial_model  489.1185  495.8223
    ## 4  penalized_spline penalized_spline_model  256.7061  271.3522

``` r
# provide models without name
compare_models(model1, model2, model3, model4)
```

    ##   model                   type       AIC       BIC
    ## 1     1       polynomial_model 1368.9096 1373.9280
    ## 2     2       polynomial_model  505.5269  508.8787
    ## 3     3       polynomial_model  489.1185  495.8223
    ## 4     4 penalized_spline_model  256.7061  271.3522

``` r
# user can provide arbitrary number of models
compare_models(model3, model4)
```

    ##   model                   type      AIC      BIC
    ## 1     1       polynomial_model 489.1185 495.8223
    ## 2     2 penalized_spline_model 256.7061 271.3522
