# Hepatitis B serological data from Russia in 1999 (aggregated)

A seroprevalence study conducted in St. Petersburg (more information in
the book)

## Usage

``` r
hbv_ru_1999
```

## Format

A data frame with 4 variables:

- age:

  Age group

- pos:

  Number of seropositive individuals

- tot:

  Total number of individuals surveyed

- gender:

  Gender of cohort (unsure what 1 and 2 means)

## Source

Mukomolov, S., L. Shliakhtenko, I. Levakova, and E. Shargorodskaya.
Viral hepatitis in Russian federation. An analytical overview. Technical
Report 213 (3), 3rd edn. St Petersburg Pasteur Institute, St Petersburg,
2000.

## Examples

``` r
library(dplyr)
hbv_ru_1999$age <- trunc(hbv_ru_1999$age / 1) * 1
hbv_ru_1999$age[hbv_ru_1999$age > 40] <- trunc(
  hbv_ru_1999$age[hbv_ru_1999$age > 40] / 5
) * 5
df <- hbv_ru_1999 %>%
  group_by(age) %>%
  summarise(pos = sum(pos), tot = sum(tot))

plot(
  df$age, df$pos / df$tot,
  cex = 0.32 * sqrt(df$tot), pch = 16, xlab = "age",
  ylab = "seroprevalence", xlim = c(0, 72)
)

```
