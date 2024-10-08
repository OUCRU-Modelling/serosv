% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weibull_model.R
\name{weibull_model}
\alias{weibull_model}
\title{The Weibull model.}
\usage{
weibull_model(t, status = NULL, pos = NULL, tot = NULL)
}
\arguments{
\item{t}{the time vector.}

\item{status}{the serostatus vector (optional if pos & tot are provided).}

\item{pos}{the positive count vector (optional if status is provided).}

\item{tot}{the total count vector (optional if status is provided).}
}
\value{
list of class weibull_model with the following items
  \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
  \item{df}{the dataframe used for fitting the model}
  \item{info}{fitted "glm" object}
  \item{sp}{seroprevalence}
  \item{foi}{force of infection}
}
\description{
Refers to section 6.1.2.
}
\examples{
df <- hcv_be_2006[order(hcv_be_2006$dur), ]
model <- weibull_model(
  t=df$dur,
  status=df$seropositive
  )
plot(model)

}
\seealso{
[stats::glm()] for more information on the fitted "glm" object
}
