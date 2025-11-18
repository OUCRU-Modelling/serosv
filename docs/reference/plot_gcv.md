# Plotting GCV values with respect to different nn-s and h-s parameters.

Refers to section 7.2.

## Usage

``` r
plot_gcv(age, pos, tot, nn_seq, h_seq, kern = "tcub", deg = 2)
```

## Arguments

- age:

  the age vector.

- pos:

  the pos vector.

- tot:

  the tot vector.#'

- nn_seq:

  Nearest neighbor sequence.

- h_seq:

  Smoothing parameter sequence.

- kern:

  Weight function, default = "tcub". Other choices are "rect", "trwt",
  "tria", "epan", "bisq" and "gauss". Choices may be restricted when
  derivatives are required; e.g. for confidence bands and some bandwidth
  selectors.

- deg:

  Degree of polynomial to use. Default: 2.

## Value

plot of gcv value

## Examples

``` r
df <- mumps_uk_1986_1987
plot_gcv(
  df$age, df$pos, df$tot,
  nn_seq = seq(0.2, 0.8, by=0.1),
  h_seq = seq(5, 25, by=1)
)

```
