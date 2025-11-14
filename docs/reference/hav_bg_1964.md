# Hepatitis A serological data from Bulgaria in 1964 (aggregated)

A cross-sectional survey conducted in 1964 in Bulgaria. Samples were
collected from schoolchildren and blood donors.

## Usage

``` r
hav_bg_1964
```

## Format

A data frame with 3 variables:

- age:

  Age group

- pos:

  Number of seropositive individuals

- tot:

  Total number of individuals surveyed

## Source

Keiding, Niels. "Age-Specific Incidence and Prevalence: A Statistical
Perspective." Journal of the Royal Statistical Society. Series A
(Statistics in Society) 154, no. 3 (1991): 371-412.
[doi:10.2307/2983150](https://doi.org/10.2307/2983150)

## Examples

``` r
# Reproduce Fig 4.1 (lower panel), p. 63
age <- hav_bg_1964$age
pos <- hav_bg_1964$pos
tot <- hav_bg_1964$tot
plot(
    age, pos / tot,
    pty = "s", cex = 0.08 * tot, pch = 16, xlab = "age",
    ylab = "seroprevalence", xlim = c(0, 86), ylim = c(0, 1)
)

```
