## Minor release
In this version I have
* add to_titer() function
* update compare_models() to allow custom model selection function
* update documentation


## R CMD check results
checking installed package size ... NOTE
installed size is  5.2Mb
sub-directories of 1Mb or more:
  doc    1.2Mb
  libs   3.1Mb
  
checking package dependencies ... NOTE
Imports includes 21 non-default packages.
Importing from so many packages makes the package vulnerable to any of
them becoming unavailable.  Move as many as possible to Suggests and
use conditionally.

0 errors | 0 warnings | 2 notes

* doc directories are for the documentation and code examples. libs contains the compiled C++ code necessary for the package.

* while we acknowledge the large number of dependencies, these imports are each directly required by 
exported functions as the package aims to support a wide range of distinct statistical model.

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

