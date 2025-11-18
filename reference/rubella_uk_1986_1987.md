# Rubella serological data from the UK in 1986 and 1987 (aggregated)

Prevalence of rubella in the UK, obtained from a large survey of
prevalence of antibodies to both mumps and rubella viruses.

## Usage

``` r
rubella_uk_1986_1987
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

Morgan-Capner P, Wright J, Miller C L, Miller E. Surveillance of
antibody to measles, mumps, and rubella by age. British Medical Journal
1988; 297 :770
[doi:10.1136/bmj.297.6651.770](https://doi.org/10.1136/bmj.297.6651.770)

## Examples

``` r
# Reproduce Fig 4.4 (middle panel), p. 67
age <- rubella_uk_1986_1987$age
pos <- rubella_uk_1986_1987$pos
tot <- rubella_uk_1986_1987$tot
plot(age, pos / tot,
  cex = 0.008 * tot, pch = 16, xlab = "age", ylab = "seroprevalence",
  xlim = c(0, 45), ylim = c(0, 1)
)

```
