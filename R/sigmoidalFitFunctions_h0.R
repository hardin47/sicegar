#' @title Sigmoidal fit function with h0 estimation
#'
#' @param dataInput A data frame or a list containing the dataframe. The data frame should be composed of at least two columns. One represents time, and the other represents intensity. The data should be normalized with the normalize data function sicegar::normalizeData() before imported into this function.
#' @param tryCounter  A counter that shows the number of times the data was fit via maximum likelihood function.
#' @param startList A vector containing the initial set of parameters that the algorithm tries for the first fit attempt for the relevant parameters. The vector composes of four elements; 'maximum', 'slopeParam', 'midPoint', and 'h0'.  Detailed explanations of those parameters can be found in vignettes. Defaults are maximum = 1, slopeParam = 1, midPoint = 0.33, and h0 = 0. The numbers are in normalized time intensity scale.
#' @param lowerBounds The lower bounds for the randomly generated start parameters.  The vector composes of four elements; 'maximum', 'slopeParam', 'midPoint', and 'h0'. Detailed explanations of those parameters can be found in vignettes. Defaults are maximum = 0.3, slopeParam = 0.01, midPoint = -0.52, and h0 = -0.1. The numbers are in normalized time intensity scale.
#' @param upperBounds The upper bounds for the randomly generated start parameters.  The vector composes of four elements; 'maximum', 'slopeParam','midPoint', and 'h0'. Detailed explanations of those parameters can be found in vignettes. Defaults are maximum = 1.5, slopeParam = 180, midPoint = 1.15, and h0 = 0.3. The numbers are in normalized time intensity scale.
#' @param min_Factor Defines the minimum step size used by the fitting algorithm. Default is 1/2^20.
#' @param n_iterations Defines maximum number of iterations used by the fitting algorithm. Default is 1000
#'
#' @description The function fits a sigmoidal curve to given data by using likelihood maximization (LM) algorithm and provides the parameters (maximum, slopeParam, midPoint, and h0) describing the double-sigmoidal fit as output. It also contains information about the goodness of fits such as AIC, BIC, residual sum of squares, and log likelihood.
#' @return Returns fitted parameters for the sigmoidal model.
#' @export
#'
#' @examples
#'time <- seq(3, 24, 0.5)
#'
#'#simulate intensity data and add noise
#'noise_parameter <- 0.1
#'intensity_noise <- stats::runif(n = length(time), min = 0, max = 1) * noise_parameter
#'intensity <- sigmoidalFitFormula_h0(time, maximum = 4, slopeParam = 1, midPoint = 8, h0 = 1)
#'intensity <- intensity + intensity_noise
#'
#'dataInput <- data.frame(intensity = intensity, time = time)
#'normalizedInput <- normalizeData(dataInput)
#'parameterVector <- sigmoidalFitFunction_h0(normalizedInput, tryCounter = 1)
#'
#'#Check the results
#'# sigmoidalFitFunction_h0() is run on the startList param values (because 'tryCounter = 1')
#'# use multipleFitFunction() for multiple random starts in order to optimize
#'if(parameterVector$isThisaFit){
#'intensityTheoretical <- sigmoidalFitFormula_h0(time,
#'                                            maximum = parameterVector$maximum_Estimate,
#'                                            slopeParam = parameterVector$slopeParam_Estimate,
#'                                            midPoint = parameterVector$midPoint_Estimate,
#'                                            h0 = parameterVector$h0_Estimate)
#'
#'comparisonData <- cbind(dataInput, intensityTheoretical)
#'
#'require(ggplot2)
#'ggplot(comparisonData) +
#'  geom_point(aes(x = time, y = intensity)) +
#'  geom_line(aes(x = time, y = intensityTheoretical)) +
#'  expand_limits(x = 0, y = 0)
#'}
#'
#'if(!parameterVector$isThisaFit){
#'   print(parameterVector)
#'}
#'
sigmoidalFitFunction_h0 <- function(dataInput,
                                 tryCounter,
                                 startList = list(maximum = 1, slopeParam = 1, midPoint = 0.33, h0 = 0), #adding bounds for h0
                                 lowerBounds = c(maximum = 0.3, slopeParam = 0.01,  midPoint = -0.52, h0 = -0.1),
                                 upperBounds = c(maximum = 1.5, slopeParam = 180,  midPoint = 1.15, h0 = 0.3),
                                 min_Factor = 1/2^20,
                                 n_iterations = 1000){

  isalist <- (is.list(dataInput) & !is.data.frame(dataInput))
  if(isalist){
    dataFrameInput <- dataInput$timeIntensityData
  }

  isadataframe <- (is.data.frame(dataInput))

  if(isadataframe){
    dataFrameInput <- dataInput
  }

  if(tryCounter == 1){
    counterDependentStartList <- startList
  }

  if(tryCounter != 1){
    randomVector <- stats::runif(length(startList), 0, 1)
    names(randomVector) <- c("maximum", "slopeParam", "midPoint", "h0") #added h0
    counterDependentStartVector <- randomVector * (upperBounds - lowerBounds) + lowerBounds
    counterDependentStartList <- as.list(counterDependentStartVector)
  }

  theFitResult <- try(minpack.lm::nlsLM(intensity ~ sicegar::sigmoidalFitFormula_h0(time, maximum, slopeParam, midPoint, h0), #added h0 argument
                                        dataFrameInput,
                                        start = counterDependentStartList,
                                        control = list(maxiter = n_iterations, minFactor = min_Factor),
                                        lower = lowerBounds,
                                        upper = upperBounds,
                                        trace = F), silent = TRUE)

  if(!inherits(theFitResult, "try-error")){

    parameterMatrix <- summary(theFitResult)$parameters
    colnames(parameterMatrix) <- c("Estimate", "Std_Error", "t_value", "Pr_t")

    parameterVector <- c(t(parameterMatrix))
    names(parameterVector) <- c("maximum_N_Estimate", "maximum_Std_Error", "maximum_t_value", "maximum_Pr_t",
                                "slopeParam_N_Estimate", "slopeParam_Std_Error", "slopeParam_t_value", "slopeParam_Pr_t",
                                "midPoint_N_Estimate", "midPoint_Std_Error", "midPoint_t_value", "midPoint_Pr_t",
                                "h0_N_Estimate",    "h0_Std_Error",    "h0_t_value",    "h0_Pr_t") #added h0

    parameterVector <- c(parameterVector,
                         residual_Sum_of_Squares = sum((as.vector(stats::resid(theFitResult)))^2)[1],
                         log_likelihood = as.vector(stats::logLik(theFitResult))[1],
                         AIC_value = as.vector(stats::AIC(theFitResult))[1],
                         BIC_value = as.vector(stats::BIC(theFitResult))[1])

    parameterList <- as.list(parameterVector)
    parameterList$isThisaFit <- TRUE
    parameterList$startVector <- counterDependentStartList

    if(isalist){
      parameterList$dataScalingParameters <- as.list(dataInput$dataScalingParameters)
    }

    parameterList$model <- as.character("sigmoidal")
    parameterList$additionalParameters <- FALSE

    parameterDf <- as.data.frame(parameterList)

    #Renormalize Parameters
    parameterDf <- sigmoidalRenormalizeParameters_h0(parameterDf, isalist) #changed to h0 function

  }

  if(inherits(theFitResult, "try-error")){

    parameterVector=rep(NA, 16) #changed from 12 to 16
    names(parameterVector) <- c("maximum_N_Estimate", "maximum_Std_Error", "maximum_t_value", "maximum_Pr_t",
                                "slopeParam_N_Estimate", "slopeParam_Std_Error", "slopeParam_t_value", "slopeParam_Pr_t",
                                "midPoint_N_Estimate", "midPoint_Std_Error", "midPoint_t_value", "midPoint_Pr_t",
                                "h0_N_Estimate",    "h0_Std_Error",    "h0_t_value",    "h0_Pr_t")  #added h0

    parameterVector <- c(parameterVector,
                         residual_Sum_of_Squares = Inf,
                         log_likelihood = NA,
                         AIC_value = NA,
                         BIC_value = NA)

    parameterList <- as.list(parameterVector)
    parameterList$isThisaFit <- FALSE
    parameterList$startVector <- counterDependentStartList

    if(isalist){
      parameterList$dataScalingParameters <- as.list(dataInput$dataScalingParameters)
    }

    parameterList$model <- "sigmoidal"

    parameterDf <- as.data.frame(parameterList)

    #Renormalize Parameters
    parameterDf <- sigmoidalRenormalizeParameters_h0(parameterDf, isalist) #changed to h0 function

  }

  return(parameterDf)
}

#**************************************
#' @title sigmoidalFitFormula_h0
#'
#' @param x  the "time" (time) column of the dataframe.
#' @param maximum the maximum intensity that the sigmoidal function can reach while time approaches infinity.
#' @param slopeParam  the slope parameter of the sigmoidal function at the steepest point.
#' @param midPoint  the x axis value of the steepest point in the function.
#' @param h0 the lower asymptote (baseline) intensity
#'
#' @description Calculates intensities for given time points (x) by using sigmoidal fit model and parameters (maximum, slopeParam, midpoint, and h0).
#' @return Returns the predicted intensities for given time points with the given sigmoidal fit parameters.
#'
#' @examples
#'
#'time <- seq(3, 24, 0.5)
#'
#'#simulate intensity data and add noise
#'noise_parameter <- 0.1
#'intensity_noise <- stats::runif(n = length(time), min = 0, max = 1) * noise_parameter
#'intensity <- sigmoidalFitFormula_h0(time, maximum = 4, slopeParam = 1, midPoint = 8, h0 = 1)
#'intensity <- intensity + intensity_noise
#'
#'dataInput <- data.frame(intensity = intensity, time = time)
#'normalizedInput <- normalizeData(dataInput)
#'parameterVector <- sigmoidalFitFunction_h0(normalizedInput, tryCounter = 1)
#'
#'#Check the results
#'# sigmoidalFitFunction_h0() is run on the startList param values (because 'tryCounter = 1')
#'# use multipleFitFunction() for multiple random starts in order to optimize
#'if(parameterVector$isThisaFit){
#'  intensityTheoretical <- sigmoidalFitFormula_h0(time,
#'                                              maximum = parameterVector$maximum_Estimate,
#'                                              slopeParam = parameterVector$slopeParam_Estimate,
#'                                              midPoint = parameterVector$midPoint_Estimate,
#'                                              h0 = parameterVector$h0_Estimate)
#'
#'  comparisonData <- cbind(dataInput, intensityTheoretical)
#'
#'  require(ggplot2)
#'  ggplot(comparisonData) +
#'    geom_point(aes(x = time, y = intensity)) +
#'    geom_line(aes(x = time, y = intensityTheoretical)) +
#'    expand_limits(x = 0, y = 0)
#'}
#'
#'if(!parameterVector$isThisaFit){
#'   print(parameterVector)
#'}
#'
#'
#'
#' @export
sigmoidalFitFormula_h0 <- function(x, maximum, slopeParam, midPoint, h0){
  y=(h0 + (maximum - h0)/(1 + exp((-slopeParam) * (x - midPoint))));
  return(y)
}
#**************************************


#**************************************
# @title sigmoidalRenormalizeParameters_h0 (This is an internal function)
# @param parametersDf it is the parameter data frame generated by sigmoidal fit function
#        includes the parameters named
#        *maximum_N_Estimate (normalized Maximum Estimate)
#        *slopeParam_N_Estimate (normalzied Slope Parameter Estimate)
#        *midPoint_N_Estimate (normalized Midpoint Estimate)
#        *h0_N_Estimate (normalized h0 Estimate)
#        *dataScalingParameters.intensityRange the range of initial unnormalized intensity. Provided if the data is normalized
#        *parameterDF$dataScalingParameters.intensityMin the minimum of unnormalized intensity. Provided if the data is normalized
#        *parameterDF$dataScalingParameters.timeRange th maximum time that the experiment reach. Provided if the data is normalized
# @param isalist defines if the input is provided as a list (i.e normalized) or as a data frame (i.e not normalized)
# @details If the fit was done in normalized data frame then the found fit parameters will belong to normalized data.
#          This function generates unnormalized counterparts of those parameters.
sigmoidalRenormalizeParameters_h0 <- function(parameterDF, isalist){

  if(isalist){
    parameterDF$maximum_Estimate <- (parameterDF$maximum_N_Estimate * parameterDF$dataScalingParameters.intensityRange) + parameterDF$dataScalingParameters.intensityMin
    parameterDF$slopeParam_Estimate <- parameterDF$slopeParam_N_Estimate / parameterDF$dataScalingParameters.timeRange
    parameterDF$midPoint_Estimate <- parameterDF$midPoint_N_Estimate * parameterDF$dataScalingParameters.timeRange

    #adding h0
    parameterDF$h0_Estimate <- (parameterDF$h0_N_Estimate * parameterDF$dataScalingParameters.intensityRange) + parameterDF$dataScalingParameters.intensityMin
  }

  if(!isalist){
    parameterDF$maximum_Estimate <- parameterDF$maximum_N_Estimate
    parameterDF$slopeParam_Estimate <- parameterDF$slopeParam_N_Estimate
    parameterDF$midPoint_Estimate <- parameterDF$midPoint_N_Estimate

    #adding h0
    parameterDF$h0_Estimate <- parameterDF$h0_N_Estimate
  }

  return(parameterDF)
}
#**************************************
