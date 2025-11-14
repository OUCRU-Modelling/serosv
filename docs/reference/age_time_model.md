# Age-time varying seroprevalence

Fit age-stratified seroprevalence across multiple time points. Also try
to monotonize age (or birth cohort) - specific seroprevalence.

## Usage

``` r
age_time_model(
  data,
  time_col = "date",
  grouping_col = "group",
  age_correct = F,
  le = 512,
  ci = 0.95,
  monotonize_method = "pava"
)
```

## Arguments

- data:

  \- input data, must have\`age\`, \`status\`, time, group columns,
  where group column determines how data is aggregated

- time_col:

  \- name of the column for time (default to \`date\`)

- grouping_col:

  \- name of the column for time (default to \`group\`)

- age_correct:

  \- a boolean, if \`TRUE\`, monotonize age-specific prevalence.
  Monotonize birth cohort-specific seroprevalence otherwise.

- le:

  \- number of bins to generate age grid, used when monotonizing data

- ci:

  \- confidence interval for smoothing

- monotonize_method:

  \- either "pava" or "scam"

## Value

a list of class time_age_model with 4 items

- out:

  a data.frame with dimension n_group x 9, where columns \`info\`,
  \`sp\`, \`foi\` store output for non-monotonized data and
  \`monotonized_info\`, \`monotonized_sp\`, \`monotonized_foi\`,
  \`monotonized_ci_mod\` store output for monotonized data

- grouping_col:

  name of the column for grouping

- age_correct:

  a boolean indicating whether the data is monotonized across age or
  cohort

- datatype:

  whether the input data is aggregated or line-listing data
