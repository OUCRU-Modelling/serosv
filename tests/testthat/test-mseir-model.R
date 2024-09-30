test_that("mseir_model returns expected results", {
  expected = list(
    a=4.000000e+01,
    m=1.804851e-35,
    s=3.727363e-04,
    e=2.881610e-06,
    i=2.069520e-06,
    r=9.996223e-01)

  # last.age,     d, lambda, sigma, ni
  #       40, 1/0.5,    0.2, 26.07, 36.5
  output <- mseir_model(
    a=seq(from=1,to=40,length=500), # age range from 0 -> 40 yo
    gamma=1/0.5, # 6 months in the maternal antibodies
    lambda=0.2,  # 5 years in the susceptible class
    sigma=26.07, # 14 days in the latent class
    nu=36.5      # 10 days in the infected class
  )$output
  actual <- as.list(tail(output, 1))

  expect_equal(actual, expected, tolerance = 0.000001)
})
