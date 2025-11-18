# A bivariate Dale model

Refers to section 12.2

## Usage

``` r
bivariate_dale_model(age, y, constant_or = F, monotonized = T, sparopt = 0)
```

## Arguments

- age:

  the age vector

- y:

  a dataframe of predictors (NN, NP, PN, PP), does not have to be in
  order

- constant_or:

  whether to set constant OR restriction or not

- monotonized:

  whether to monotonize seroprevalence or not

- sparopt:

  smoothing parameter for VGAM model

## Examples

``` r
data <- rubella_mumps_uk
model <- bivariate_dale_model(age = data$age,
        y = data[, c("NN", "NP", "PN", "PP")], monotonized=TRUE)
plot(model, y1 = "Rubella", y2 = "Mumps", plot_type = "sp")

```
