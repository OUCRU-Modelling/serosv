## R CMD check results
Examples with CPU (user + system) or elapsed time > 10s
                             user system elapsed
hierarchical_bayesian_model 29.17   0.03   29.43
sir_subpops_model           12.00   0.66   12.65

checking installed package size ... NOTE
installed size is  5.2Mb
sub-directories of 1Mb or more:
  doc    1.2Mb
  libs   3.1Mb

0 errors | 0 warnings | 2 notes

* Function hierarchical_bayesian_model and sir_subpops_model are supposed to have a long run time due to sampling process and iterative nature of the function
* doc directories are for the documentation and code examples. libs contains the compiled C++ code necessary for the package.

* This is a new release.

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

