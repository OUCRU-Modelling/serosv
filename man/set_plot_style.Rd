% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{set_plot_style}
\alias{set_plot_style}
\title{Helper to adjust styling of a plot}
\usage{
set_plot_style(
  sero = "blueviolet",
  ci = "royalblue1",
  foi = "#fc0328",
  sero_line = "solid",
  foi_line = "dashed",
  xlabel = "Age"
)
}
\arguments{
\item{sero}{- color for seroprevalence line}

\item{ci}{- color for confidence interval}

\item{foi}{- color for force of infection line}

\item{sero_line}{- linetype for seroprevalence line}

\item{foi_line}{- linetype for force of infection line}

\item{xlabel}{- x label}
}
\value{
list of updated aesthetic values
}
\description{
Helper to adjust styling of a plot
}
