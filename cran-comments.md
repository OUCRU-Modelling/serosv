## R CMD check results

* checking R code for possible problems ... [24s] NOTE
compute_ci: no visible binding for global variable '.'
compute_ci: no visible binding for global variable 'fit'
compute_ci: no visible binding for global variable 'se.fit'
compute_ci.fp_model: no visible binding for global variable '.'
compute_ci.fp_model: no visible binding for global variable 'fit'
compute_ci.fp_model: no visible binding for global variable 'se.fit'
compute_ci.penalized_spline_model: no visible binding for global
  variable '.'
compute_ci.penalized_spline_model: no visible binding for global
  variable 'fit'
compute_ci.penalized_spline_model: no visible binding for global
  variable 'se.fit'
compute_ci.weibull_model: no visible binding for global variable '.'
compute_ci.weibull_model: no visible binding for global variable 'fit'
compute_ci.weibull_model: no visible binding for global variable
  'se.fit'
plot.estimate_from_mixture: no visible binding for global variable
  'pos'
plot.estimate_from_mixture: no visible binding for global variable
  'tot'
plot.estimate_from_mixture: no visible binding for global variable
  'foi_x'
plot.estimate_from_mixture: no visible binding for global variable
  'foi'
  
Examples with CPU (user + system) or elapsed time > 10s
                             user system elapsed
hierarchical_bayesian_model 29.17   0.03   29.43
sir_subpops_model           12.00   0.66   12.65

0 errors | 0 warnings | 2 note

* The variables said to have no visible binding are variables of input dataframe accessed via tidyverse syntax. As the input must be the specified object, the dataframe should always have these variables.

* Function hierarchical_bayesian_model and sir_subpops_model are supposed to have a long run time due to sampling process and iterative nature of the function

* This is a new release.

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

