% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correct_prevalence.R
\name{correct_prevalence}
\alias{correct_prevalence}
\title{Estimate the true sero prevalence using Bayesian estimation}
\usage{
correct_prevalence(
  data,
  bayesian = TRUE,
  init_se = 0.95,
  init_sp = 0.8,
  study_size_se = 1000,
  study_size_sp = 1000,
  chains = 1,
  warmup = 1000,
  iter = 2000
)
}
\arguments{
\item{data}{the input data frame, must either have `age`, `pos`, `tot` columns (for aggregated data) OR `age`, `status` for (linelisting data)}

\item{bayesian}{whether to adjust sero-prevalence using the Bayesian or frequentist approach. If set to `TRUE`, true sero-prevalence is estimated using MCMC.}

\item{init_se}{sensitivity of the serological test}

\item{init_sp}{specificity of the serological test}

\item{study_size_se}{(applicable when `bayesian=TRUE`) study size for sensitivity validation study (i.e., number of confirmed infected patients in the study)}

\item{study_size_sp}{(applicable when `bayesian=TRUE`) study size for specificity validation study (i.e., number of confirmed non-infected patients in the study)}

\item{chains}{(applicable when `bayesian=TRUE`) number of Markov chains}

\item{warmup}{(applicable when `bayesian=TRUE`) number of warm up runs}

\item{iter}{(applicable when `bayesian=TRUE`) number of iterations}
}
\value{
a list of 2 items
  \item{info}{estimated parameters}
  \item{corrected_sero}{data.frame containing age, the corresponding estimated seroprevalance, adjusted tot and pos}
}
\description{
Estimate the true sero prevalence using Bayesian estimation
}
\examples{
data <- rubella_uk_1986_1987
correct_prevalence(data)
}
