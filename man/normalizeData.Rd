% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalizationFunction.R
\name{normalizeData}
\alias{normalizeData}
\title{Normalization of given data}
\usage{
normalizeData(dataInput, dataInputName = NA)
}
\arguments{
\item{dataInput}{A data frame or a list contatining the dataframe. The data frame should be composed of at least two columns. One represents time, and the other represents intensity.}

\item{dataInputName}{experiment name (Default is 'NA').}
}
\value{
Function returns a new data frame, scaling factors and scaling constants that connects initial data frame to new one. The new data frame includes 2 columns one is for normalized time and the other is for noralized intensity. The whole time is distributed between 0 and 1 and similarly the whole intensity is distributed between 0 and 1. The time and intensity constants and scaling factors are the parameters to transform data from unnormalized data frame to normalized data frame.
}
\description{
Maps the given time-intensity data into a rescaled data frame where time is scaled in a way that maximum time point is one and intensity is distributed between [0,1].
}
\examples{
# generateRandomData
time <- seq(1, 48, 0.5)
intensity <- runif(length(time), 3.0, 7.5)
dataInput <- data.frame(time, intensity)

# Normalize Data
dataOutput <- normalizeData(dataInput, dataInputName="sample001")

}
