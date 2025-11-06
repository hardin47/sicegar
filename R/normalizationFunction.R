#' @title Normalization of given data
#'
#' @param dataInput  A dataframe or a list containing the dataframe. The data frame should be composed of at least two columns. One represents time, and the other represents intensity.
#' @param dataInputName experiment name (Default is 'NA').
#' @return Function returns a new dataframe, scaling factors and scaling constants that connects the initial data frame to the new one. The new data frame includes 2 columns: normalized time and normalized intensity. The time and intensity constants and scaling factors are the parameters to transform data from the unnormalized dataframe to normalized data frame.
#' @description Maps the given time-intensity data into a rescaled dataframe where time is scaled to between 0 and 1, and intensity is scaled to be between 0 and 1.
#'
#' @export
#'
#' @examples
#' # runif() is used here for consistency with previous versions of the sicegar package. However,
#' # rnorm() will generate symmetric errors, which will produce less biased numerical estimates of the parameters.
#' # We recommend errors generated with rnorm() for any simulation studies on sicegar.
#' # generateRandomData
#' time <- seq(3, 48, 0.5)
#' intensity <- runif(length(time), 3.0, 7.5)
#' dataInput <- data.frame(time, intensity)
#'
#' # Normalize Data
#' dataOutput <- normalizeData(dataInput, dataInputName="sample001")
#'
normalizeData <-function(dataInput, dataInputName = NA)
  {
    dataInputCheckVariable <- dataCheck(dataInput)

    #timeMin <- min(dataInput$time)
    timeData <- dataInput$time
    timeRange <- max(timeData,na.rm = T)
    timeData <- timeData / timeRange

    intensityMin <- min(dataInput$intensity,na.rm = T)
    intensityMax <- max(dataInput$intensity,na.rm = T)
    intensityData <- dataInput$intensity - intensityMin
    intensityRange <- max(intensityData,na.rm = T)
    intensityData <- intensityData / intensityRange

    dataOutput <- data.frame(time = timeData, intensity = intensityData)
    return(list(timeIntensityData = dataOutput,
                dataScalingParameters = c(timeRange = timeRange,
                                          intensityMin = intensityMin,
                                          intensityMax = intensityMax,
                                          intensityRange = intensityRange),
                dataInputName = dataInputName))
  }


#' @title Unnormalization of given data
#'
#' @param dataInput a list file composed of two parts
#' First part is the data that will be unnormalized, which is a dataframe composed of two columns. One is for time and the other is for intensity
#' Second part is the scaling parameters of the data which is a vector that has three components. The first is related with time and second two are related with intensity. The second value represents the min value of the intensity set. First and third values represent the difference between max and min value in the relevant column.
#' @return Returns a dataframe, scaling factors, and scaling constants for time and intensity. The other data frame includes 2 columns: normalized time and normalized intensity. The time and intensity constants and scaling factors are the parameters to transform data from given set to scaled set.
#' @description Maps the given time-intensity data into a rescaled frame where time is between zero and one intensity is also between zero and one.
#'
#' @export
#'
#' @examples
#' # runif() is used here for consistency with previous versions of the sicegar package. However,
#' # rnorm() will generate symmetric errors, which will produce less biased numerical estimates of the parameters.
#' # We recommend errors generated with rnorm() for any simulation studies on sicegar.
#' # generateRandomData
#' time <- seq(3, 48, 0.5)
#' intensity <- runif(length(time), 3.0, 7.5)
#' dataInput <- data.frame(time, intensity)
#' # Normalize Data
#' dataOutput <- normalizeData(dataInput)
#' dataInput2 <- dataOutput
#' # Un Normalize it
#' dataOutput2 <- unnormalizeData(dataInput2)
#'
unnormalizeData <-
  function(dataInput)
  {
    dataInputCheckVariable <- dataCheck(dataInput)
    time <- dataInput$dataScalingParameters[["timeRange"]] * dataInput$timeIntensityData[["time"]]
    intensity <- dataInput$dataScalingParameters[["intensityMin"]] +
      dataInput$dataScalingParameters[["intensityRange"]] * dataInput$timeIntensityData[["intensity"]]

    dataOutput <- list(timeIntensityData = data.frame(time, intensity))
    return(dataOutput)
  }

