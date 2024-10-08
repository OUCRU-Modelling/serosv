% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_ci.R
\name{compute_ci.fp_model}
\alias{compute_ci.fp_model}
\title{Compute confidence interval for fractional polynomial model}
\usage{
compute_ci.fp_model(x, ci = 0.95, le = 100, ...)
}
\arguments{
\item{x}{- serosv models}

\item{ci}{- confidence interval}

\item{le}{- number of data for computing confidence interval}

\item{...}{- arbitrary argument}
}
\value{
confidence interval dataframe with 4 variables, x and y for the fitted values and ymin and ymax for the confidence interval
}
\description{
Compute confidence interval for fractional polynomial model
}
