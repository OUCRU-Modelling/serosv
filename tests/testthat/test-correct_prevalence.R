test_that("test methods for correcting prevalence", {
  data <- rubella_uk_1986_1987

  # test bayesian method
  expect_no_error(
    out_bayesian <- correct_prevalence(
      data,
      warmup = 500,
      iter = 2000,
      init_se = 0.9,
      init_sp = 0.8,
      study_size_se = 1000,
      study_size_sp = 3000
    )
  )

  # test frequentist method
  expect_no_error(
    out_freq <- correct_prevalence(data, bayesian = FALSE, init_se=0.9, init_sp = 0.8)
  )


  # test plot functions for both approaches
  expect_no_error(plot_corrected_prev(out_freq))
  expect_no_error(plot_corrected_prev(out_bayesian))

  # making sure that facet settings also works
  expect_no_error(plot_corrected_prev(out_bayesian, out_freq, facet = TRUE))
  expect_no_error(plot_corrected_prev(out_bayesian, out_freq, facet = FALSE))
})
