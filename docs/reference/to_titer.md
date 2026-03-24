# Convert assay readings to titers

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

  - A string naming a built-in model (currently supported: \`"4PL"\`),
    or

  - A named list with two functions: \`\$mod\` for curve fitting and
    \`\$quantify_ci\` for titer estimation with confidence intervals.

- positive_threshold:

  if not NULL, processed_data will have the serostatus labeled

- ci:

  confidence interval for the titer estimates (default is .95 i.e., 95

  negative_controlif TRUE, output tibble will include the result for
  negative controls

a data.frame with 8 columns to_titer() converts raw assay readings
(e.g., OD, fluorescence intensity) to titer by fitting a calibrating
model
