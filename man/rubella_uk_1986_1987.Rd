% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{rubella_uk_1986_1987}
\alias{rubella_uk_1986_1987}
\title{Rubella serological data from the UK in 1986 and 1987 (aggregated)}
\format{
A data frame with 3 variables:
\describe{
 \item{age}{Age group}
 \item{pos}{Number of seropositive individuals}
 \item{tot}{Total number of individuals surveyed}
}
}
\source{
Morgan-Capner P, Wright J, Miller C L, Miller E. Surveillance of
antibody to measles, mumps, and rubella by age. British Medical Journal
1988; 297 :770 \doi{doi:10.1136/bmj.297.6651.770}
}
\usage{
rubella_uk_1986_1987
}
\description{
Prevalence of rubella in the UK, obtained from a large survey of prevalence
of antibodies to both mumps and rubella viruses.
}
\examples{
# Reproduce Fig 4.4 (middle panel), p. 67
age <- rubella_uk_1986_1987$age
pos <- rubella_uk_1986_1987$pos
tot <- rubella_uk_1986_1987$tot
plot(age, pos / tot,
  cex = 0.008 * tot, pch = 16, xlab = "age", ylab = "seroprevalence",
  xlim = c(0, 45), ylim = c(0, 1)
)

}
\keyword{datasets}
