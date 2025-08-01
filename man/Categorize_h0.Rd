% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/categorize_h0.R
\name{Categorize_h0}
\alias{Categorize_h0}
\title{Decide between Sigmoidal and Double-Sigmoidal Fits}
\usage{
Categorize_h0(
  parameterVectorSigmoidal,
  parameterVectorDoubleSigmoidal,
  threshold_intensity_range = 0.1,
  threshold_minimum_for_intensity_maximum = 0.3,
  threshold_bonus_sigmoidal_AIC = 0,
  threshold_sm_tmax_IntensityRatio = 0.75,
  threshold_dsm_tmax_IntensityRatio = 0.75,
  threshold_AIC = -10,
  threshold_t0_max_int = 0.05,
  showDetails = FALSE
)
}
\arguments{
\item{parameterVectorSigmoidal}{A named list or data frame of sigmoidal fit parameters,
produced by \code{\link{sigFit}} and
\code{\link{parameterCalculation_h0}}.}

\item{parameterVectorDoubleSigmoidal}{A named list or data.frame of double-sigmoidal fit
parameters, produced by \code{\link{sigFit2}} and
\code{\link{parameterCalculation_h0}}.}

\item{threshold_intensity_range}{Numeric; minimum required intensity range to consider
any signal (default 0.1).}

\item{threshold_minimum_for_intensity_maximum}{Numeric; minimum maximum intensity
required to consider a real signal (default 0.3).}

\item{threshold_bonus_sigmoidal_AIC}{Numeric; AIC bonus given to the sigmoidal model
when comparing with the double-sigmoidal (default 0).}

\item{threshold_sm_tmax_IntensityRatio}{Numeric; minimum fraction of maximum
intensity reached by sigmoidal at final time to allow sigmoidal (default 0.85).}

\item{threshold_dsm_tmax_IntensityRatio}{Numeric; maximum fraction of maximum
intensity reached by double-sigmoidal at final time to allow double-sigmoidal (default 0.75).}

\item{threshold_AIC}{Numeric; maximum AIC value to allow either model (default –10).}

\item{threshold_t0_max_int}{Numeric; maximum allowed starting intensity ratio
at time zero (default 0.05).}

\item{showDetails}{Logical; if TRUE, prints the internal decision list structure
for debugging (default FALSE).}
}
\value{
A named list with components:
  \item{decisionList}{All intermediate TRUE/FALSE tests and thresholds.}
  \item{decision}{One of “no_signal”, “sigmoidal”, “double_sigmoidal” or “ambiguous”.}
  \item{decisonSteps}{Which numbered tests were applied (as a single string).}
}
\description{
Runs a series of logical and AIC‐based tests to choose the best model
(“sigmoidal”, “double_sigmoidal”, “no_signal” or “ambiguous”) given two
fitted parameter vectors (one sigmoidal, one double-sigmoidal).
}
