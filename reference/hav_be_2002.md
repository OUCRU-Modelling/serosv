# Hepatitis A serological data from Belgium in 2002 (line listing)

A subset of the serological dataset of Varicella-Zoster Virus (VZV) and
Parvovirus B19 in Belgium where only individuals living in Flanders were
selected

## Usage

``` r
hav_be_2002
```

## Format

A data frame with 2 variables:

- age:

  Age of individual

- seropositive:

  If the individual is seropositive or not

## Source

Thiry, N., Beutels, P., Shkedy, Z. et al. The seroepidemiology of
primary varicella-zoster virus infection in Flanders (Belgium). Eur J
Pediatr 161, 588-593 (2002).
[doi:10.1007/s00431-002-1053-2](https://doi.org/10.1007/s00431-002-1053-2)

## Examples

``` r
library(dplyr)
#> Warning: package ‘dplyr’ was built under R version 4.3.1
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
df <- hav_be_2002 %>%
  group_by(age) %>%
  summarise(pos = sum(seropositive), tot = n())

with(
  df,
  plot(
    age, pos / tot,
    pty = "s", cex = 0.34 * sqrt(tot), pch = 16, xlab = "age",
    ylab = "seroprevalence", xlim = c(0, 86), ylim = c(0, 1)
  )
)

```
