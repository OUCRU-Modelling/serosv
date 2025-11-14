# Preprocess data

Preprocess data

## Usage

``` r
standardize_data(
  df,
  plate_id_col = "PLATE_ID",
  id_col = "SAMPLE_ID",
  result_col = "RESULT",
  dilution_fct_col = "DILUTION_FACTORS",
  antitoxin_label = "Anti_toxin",
  negative_col = "^NEGATIVE_*"
)
```

## Arguments

- df:

  \- data.frame with columns for plate id, sample id, result, dilution
  factor, and (optionally) negative controls

- plate_id_col:

  \- name of the column storing plates id

- id_col:

  \- name of the column storing sample id

- result_col:

  \- name of the column storing result

- dilution_fct_col:

  \- name of the column storing dilution factors

- antitoxin_label:

  \- how antitoxin is label in the sample id column

- negative_col:

  \- regex for columns for negative controls, assumed to be a label
  followed by the dilution factor (e.g. NEGATIVE_50, NEGATIVE_100)

## Value

a standardized data.frame that can be passed to \`to_titer()\`
