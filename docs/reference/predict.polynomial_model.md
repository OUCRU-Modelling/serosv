# Prediction for serosv polynomial model

A wrapper of predict.glm for direct prediction from polynomial_model
object

## Usage

``` r
# S3 method for class 'polynomial_model'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  serosv models

- newdata:

  data.frame with age column to generate prediction

- ...:

  arbitrary argument

## Value

prediction output

## See also

\[stats::predict.glm()\] for more information on the predict function
