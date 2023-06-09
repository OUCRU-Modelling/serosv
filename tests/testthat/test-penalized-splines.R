library(gam)
library(mgcv)
source("../../R/pspline.R")

test_that("tps() returns expected results", {
  expected_knots <- c(8.071233, 14.041096, 19.375342) # page 123

  df <- vzv_be_2001_2003
  filter <- (df$age>0.5)&(df$age<40)&(!is.na(df$age))&!is.na(df$parvores)
  dff <- df[filter, ]
  dfo <- dff[order(dff$age), ]
  model <- ps(
    x=dfo$age,
    y=dfo$parvores,
    k=3, deg=2
  )
  actual_knots <- model$info$info$pen$knots

  expect_equal(actual_knots, actual_knots)
})

test_that("ss(link=logit) returns expected results", {
  expected_bic <- c(3480.80) # table 8.1

  df <- parvob19_be_2001_2003
  filter <- (df$age>0.5)&(df$age<76)&(!is.na(df$age))&!is.na(df$seropositive)
  dff <- df[filter, ]
  dfo <- dff[order(dff$age), ]

  model <- ss(
    x=dfo$age,
    y=dfo$seropositive,
    d=4,
    link="logit"
  )
  plot(model, round.by=0)

  model$df$t

  actual_bic <- unname(model$info$deviance+log(length(model$info$y))*(model$info$nl.df+2))
  actual_bic

  expect_equal(actual_bic, expected_bic, tolerance=0.01)
})

test_that("ss(link=cloglog) returns expected results", {
  expected_bic <- c(3488.43) # table 8.1

  df <- vzv_be_2001_2003
  filter <- (df$age>0.5)&(df$age<76)&(!is.na(df$age))&!is.na(df$parvores)
  dff <- df[filter, ]
  dfo <- dff[order(dff$age), ]
  model <- ss(
    x=dfo$age,
    y=dfo$parvores
  )
  actual_bic <- model$info$deviance+log(length(model$info$y))*(model$info$nl.df+2)

  expect_equal(actual_bic, expected_bic, tolerance=0.01)
})

test_that("bss(link=logit) returns expected results", {
  expected_bic <- c(3472.70) # table 8.1

  df <- vzv_be_2001_2003
  filter <- (df$age>0.5)&(df$age<76)&(!is.na(df$age))&!is.na(df$parvores)
  dff <- df[filter, ]
  dfo <- dff[order(dff$age), ]
  model <- bss(
    x=dfo$age,
    y=dfo$parvores
  )
  actual_bic <- model$info$dev+log((dim(model$info$summary.predicted)[1]))*(model$info$eff.df)

  expect_equal(actual_bic, expected_bic, tolerance=0.01)
})
