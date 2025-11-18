# Generate the (NN, NP, PN, PP) matrix from line-listing data

Generate the (NN, NP, PN, PP) matrix from line-listing data

## Usage

``` r
generate_quad_matrix(data, y1, y2, age, discrete_age = F)
```

## Arguments

- data:

  \- dataframe being transformed

- y1:

  \- serostatus for y1

- y2:

  \- serostatus for y2

- age:

  \- individual's age

- discrete_age:

  \- round age so it is discrete

## Value

transformed matrix

## Examples

``` r
data <- vzv_parvo_be
data <- generate_quad_matrix(data, parvo_res, vzv_res, age)
```
