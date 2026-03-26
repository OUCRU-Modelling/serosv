# Convert assay readings to titers

to_titer() converts raw assay readings (e.g., OD, fluorescence
intensity) to titer by fitting a calibrating model

## Usage

``` r
to_titer(
  df,
  model = "4PL",
  positive_threshold = NULL,
  ci = 0.95,
  negative_control = TRUE
)
```

## Arguments

- df:

  a standardized data.frame returned by\`standardize_data()\`

- model:

  either:

  1.  A string naming a built-in model (currently supported: `"4PL"`),
      or

  2.  A named list with two functions: `$mod` for curve fitting and
      `$quantify_ci` for titer estimation with confidence intervals.

- positive_threshold:

  if not NULL, processed_data will have the serostatus labeled

- ci:

  confidence interval for the titer estimates (default is .95 i.e., 95%
  CI)

- negative_control:

  if TRUE, output tibble will include the result for negative controls

## Value

a data.frame with 8 columns

- plate_id:

  id of the plate

- data:

  list of \`data.frame\`s containing the raw sample results from each
  plate

- antitoxin_df:

  list of \`data.frame\`s containing the raw results for antitoxins from
  each plate

- standard_curve_func:

  list of functions mapping from assay reading to titer for each plate

- std_crv_midpoint:

  midpoint of the standard curve, for qualitative analysis

- processed_data:

  list of \`tibble\`s containing samples with titer estimates (lower,
  median, upper)

- negative_control:

  list of \`tibble\`s containing negative control check results (if
  \`negative_control=TRUE\`)
