% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigmoidalFitFunctions_h0.R
\name{sigmoidalFitFormula_h0}
\alias{sigmoidalFitFormula_h0}
\title{Sigmoidal Fit Formula with Baseline}
\usage{
sigmoidalFitFormula_h0(x, maximum, slopeParam, midPoint, h0)
}
\arguments{
\item{x}{Numeric vector of time points.}

\item{maximum}{Numeric; the upper asymptote of the sigmoid.}

\item{slopeParam}{Numeric; the slope parameter (rate of rise).}

\item{midPoint}{Numeric; the time at which the sigmoid reaches halfway between \code{h0} and \code{maximum}.}

\item{h0}{Numeric; the lower asymptote (baseline) intensity.}
}
\value{
Numeric vector of model‐predicted intensities.
}
\description{
Compute the sigmoidal curve value at given time points,
  allowing a nonzero baseline (\code{h0}).
}
\examples{
time <- seq(0, 100, length.out = 50)
y <- sigmoidalFitFormula_h0(time, maximum = 10, slopeParam = 0.1, midPoint = 50, h0 = 2)
}
