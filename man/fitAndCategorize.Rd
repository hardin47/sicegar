% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainFunctions.R
\name{fitAndCategorize}
\alias{fitAndCategorize}
\title{Fit and categorize.}
\usage{
fitAndCategorize(
  dataInput,
  dataInputName = NA,
  n_runs_min_sm = 20,
  n_runs_max_sm = 500,
  n_runs_min_dsm = 20,
  n_runs_max_dsm = 500,
  showDetails = FALSE,
  startList_sm = list(h0 = 0, maximum = 1, slopeParam = 1, midPoint = 0.33),
  lowerBounds_sm = c(h0 = 0, maximum = 0.3, slopeParam = 0.01, midPoint = -0.52),
  upperBounds_sm = c(h0 = 0.3, maximum = 1.5, slopeParam = 180, midPoint = 1.15),
  min_Factor_sm = 1/2^20,
  n_iterations_sm = 1000,
  startList_dsm = list(finalAsymptoteIntensityRatio = 0, maximum = 1, slope1Param = 1,
    midPoint1Param = 0.33, slope2Param = 1, midPointDistanceParam = 0.29, h0 = 0),
  lowerBounds_dsm = c(finalAsymptoteIntensityRatio = 0, maximum = 0.3, slope1Param =
    0.01, midPoint1Param = -0.52, slope2Param = 0.01, midPointDistanceParam = 0.04, h0 =
    0),
  upperBounds_dsm = c(finalAsymptoteIntensityRatio = 1, maximum = 1.5, slope1Param = 180,
    midPoint1Param = 1.15, slope2Param = 180, midPointDistanceParam = 0.63, h0 = 0),
  min_Factor_dsm = 1/2^20,
  n_iterations_dsm = 1000,
  threshold_intensity_range = 0.1,
  threshold_minimum_for_intensity_maximum = 0.3,
  threshold_bonus_sigmoidal_AIC = 0,
  threshold_sm_tmax_IntensityRatio = 0.75,
  threshold_dsm_tmax_IntensityRatio = 0.75,
  threshold_AIC = -10,
  threshold_t0_max_int = 0.05,
  stepSize = 1e-05,
  use_h0 = FALSE,
  ...
)
}
\arguments{
\item{dataInput}{Un_normalized input data that will be fitted transferred into related functions}

\item{dataInputName}{Name of data set (Default is 'NA').}

\item{n_runs_min_sm}{This number indicates the lower limit of the successful fitting attempts for sigmoidal model. It should be smaller than the upper limit of the fitting attempts (n_runs_max_sm). Default is 20}

\item{n_runs_max_sm}{This number indicates the upper limit of the fitting attempts for sigmoidal model. Default is 500}

\item{n_runs_min_dsm}{This number indicates the lower limit of the successful fitting attempts for double sigmoidal model. It should be smaller than the upper limit of the fitting attempts (n_runs_max_dsm). Default is 20}

\item{n_runs_max_dsm}{This number indicates the upper limit of the fitting attempts for sigmoidal model for double sigmoidal model. Default is 500}

\item{showDetails}{Logical if TRUE prints details of intermediate steps of individual fits (Default is FALSE).}

\item{startList_sm}{The initial set of parameters vector that sigmoidal fit algorithm tries for the first fit attempt for the relevant parameters. The vector composes of three elements; 'maximum', 'slopeParam' and, 'midPoint'.  Detailed explanations of those parameters can be found in vignettes. Defaults are maximum = 1, slopeParam = 1 and, midPoint = 0.33. The numbers are in normalized time intensity scale.}

\item{lowerBounds_sm}{The lower bounds for the randomly generated start parameters for the sigmoidal fit. The vector composes of three elements; 'maximum', 'slopeParam' and, 'midPoint'. Detailed explanations of those parameters can be found in vignettes. Defaults are maximum = 0.3, slopeParam = 0.01, and midPoint = -0.52. The numbers are in normalized time intensity scale.}

\item{upperBounds_sm}{The upper bounds for the randomly generated start parameters for the sigmoidal fit. The vector composes of three elements; 'maximum', 'slopeParam' and, 'midPoint'. Detailed explanations of those parameters can be found in vignettes. Defaults are maximum = 1.5, slopeParam = 180,  midPoint = 1.15. The numbers are in normalized time intensity scale.}

\item{min_Factor_sm}{Defines Defines the minimum step size used by the sigmoidal fit algorithm. Default is 1/2^20.}

\item{n_iterations_sm}{Defines maximum number of iterations used by the sigmoidal fit algorithm. Default is 1000}

\item{startList_dsm}{The initial set of parameters vector that double sigmoidal fit algorithm tries for the first fit attempt for the relevant parameters. The vector composes of six elements; 'finalAsymptoteIntensityRatio', 'maximum', 'slope1Param', 'midPoint1Param' , 'slope2Param', and 'midPointDistanceParam'. Detailed explanations of those parameters can be found in vignettes. Defaults are  finalAsymptoteIntensityRatio = 0, maximum = 1, slope1Param = 1, midPoint1Param = 0.33, slope2Param = 1, and midPointDistanceParam=0.29. The numbers are in normalized time intensity scale.}

\item{lowerBounds_dsm}{The lower bounds for the randomly generated start parameters for double sigmoidal fit.  The vector composes of six elements; 'finalAsymptoteIntensityRatio', 'maximum', 'slope1Param', 'midPoint1Param' , 'slope2Param', and 'midPointDistanceParam'. Detailed explanations of those parameters can be found in vignettes. Defaults are finalAsymptoteIntensityRatio = 0, maximum = 0.3, slope1Param = .01, midPoint1Param = -0.52, slope2Param = .01, and midPointDistanceParam = 0.04. The numbers are in normalized time intensity scale.}

\item{upperBounds_dsm}{The upper bounds for the randomly generated start parameters for double sigmoidal fit.  The vector composes of six elements; 'finalAsymptoteIntensityRatio', 'maximum', 'slope1Param', 'midPoint1Param' , 'slope2Param', and 'midPointDistanceParam'. Detailed explanations of those parameters can be found in vignettes. Defaults are finalAsymptoteIntensityRatio = 1, maximum = 1.5, slope1Param = 180, midPoint1Param = 1.15, slope2Param = 180, and midPointDistanceParam = 0.63. The numbers are in normalized time intensity scale.}

\item{min_Factor_dsm}{Defines the minimum step size used by the double sigmoidal fit algorithm. Default is 1/2^20.}

\item{n_iterations_dsm}{Define maximum number of iterations used by the double sigmoidal fit algorithm. Default is 1000}

\item{threshold_intensity_range}{Minimum for intensity range, i.e. it is the lower limit for the allowed difference between the maximum and minimum of the intensities (Default is 0.1, and the values are based on actual, not the rescaled data.).}

\item{threshold_minimum_for_intensity_maximum}{Minimum allowed value for intensity maximum. (Default is 0.3, and the values are based on actual, not the rescaled data.).}

\item{threshold_bonus_sigmoidal_AIC}{Bonus AIC points for sigmoidal fit. Negative values help the sigmoidal model to win. Only helps in competition between sigmoidal and double sigmoidal fit at decision step "9", i.e. if none of the models fail in any of the tests and stay as a candidate until the last step (Default is 0).}

\item{threshold_sm_tmax_IntensityRatio}{The threshold for the minimum intensity ratio between the last observed time points intensity and theoretical maximum intensity of the sigmoidal curve. If the value is below the threshold, then the data can not be represented with the sigmoidal model. (Default is 0.85)}

\item{threshold_dsm_tmax_IntensityRatio}{The threshold for the minimum intensity ratio between the last observed time points intensity and maximum intensity of the double sigmoidal curve.  If the value is above the threshold, then the data can not be represented with the double sigmoidal model. (Default is 0.75)}

\item{threshold_AIC}{Maximum AIC values in order to have a meaningful fit (Default is -10).}

\item{threshold_t0_max_int}{Maximum allowed intensity of the fitted curve at time is equal to zero (t=0). (Default is 0.05, and the values are based on actual, not the rescaled data.).}

\item{stepSize}{Step size used by the fitting algorithm. Smaller numbers gave more accurate results than larger numbers, and larger numbers gave the results faster than small numbers. The default value is 0.00001.}

\item{use_h0}{Boolean which decides whether to fix the lower asymptote of h0 at 0 (FALSE, default) or to freely estimate h0 (TRUE)}

\item{...}{All other arguments that model functions ("sigmoidalFitFunction" and, "doublesigmoidalFitFunction") may need.}
}
\value{
Returns the parameters related with the curve fitted to the input data.
}
\description{
Fits the sigmoidal and double-sigmoidal models to the data and then categorizes the data according to which model fits best.
}
\examples{
# Example 1 (For h0 fixed at zero)
time <- seq(1, 24, 0.5)

#simulate intensity data and add noise
noise_parameter <- 0.2
intensity_noise <- stats::runif(n = length(time), min = 0, max = 1) * noise_parameter
intensity <- sicegar::doublesigmoidalFitFormula(time,
                                               finalAsymptoteIntensityRatio = .3,
                                               maximum = 4,
                                               slope1Param = 1,
                                               midPoint1Param = 7,
                                               slope2Param = 1,
                                               midPointDistanceParam = 8)
intensity <- intensity + intensity_noise

dataInput <- data.frame(intensity = intensity, time = time)

fitObj1 <- sicegar::fitAndCategorize(dataInput = dataInput)

# Example 2 (For free h0)
# using the same data as above `dataInput`

fitObj2 <- sicegar::fitAndCategorize(dataInput = dataInput, use_h0 = TRUE)

}
