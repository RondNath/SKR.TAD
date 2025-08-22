library(devtools) #devtools_2.4.5
devtools::install_gitlab(repo = "urep/dev_utils/r_utils/r4urep",
                         host = "https://forge.inrae.fr")

## Step 2. Function: SKR ANALYSIS ----

#' @title Launch the analysis of the TADs with the SKR framework
#' @description Launch the SKR analysis of the TADs, and generate output dataset
#' @param weights the dataframe of weights, one row correspond to a series of observation
#' @param weightsFactor the dataframe which contains the different factor linked to the weights
#' @param dataToTreat a vector of the data linked to the different factor
#' @param randomizationNumber the number of random abundance matrix to generate
#' @param aggregationFactorName vector of factor name for the generation of random matrix
#' @param statisticsFactorName vector of factor name for the computation of statistics for each generated matrix
#' @param seed the seed of the pseudo random number generator
#' @param abundanceDataFrameRDS the path and name of the RDS file to load/save the dataframe which contains the observed data and the generated matrix
#' @param weightedMomentsDataFrameRDS the path and name of the RDS file to load/save the dataframe which contains the calculated moments
#' @param statPerObsDataFrameRDS the path and name of the RDS file to load/save the dataframe which contains the statistics for each observed row regarding the random ones
#' @param statPerRandDataFrameRDS the path and name of the RDS file to load/save the dataframe which contains the statistics for each random matrix generated
#' @param statSKRparam the path and name of the RDS file to load/save the dataframe which contains the SKR statistics
#' @param regenerateAbundanceDataFrame boolean to specify if the abundance dataframe is computed again
#' @param regenerateWeightedMomentsDataFrame boolean to specify if the weighted moments dataframe is computed again
#' @param regenerateStatPerObsDataFrame boolean to specify if the statistics per observation dataframe is computed again
#' @param regenerateStatPerRandDataFrame boolean to specify if the statistics per random matrix dataframe is computed again
#' @param significanceThreshold the significance threshold to consider that the observed value is in the randomized value
#' @param doParallel Indicates if we use parallelism to construct the random matrix
#' @param slope_speTADs slope of a specific SKR used as a baseline (default: slope_speTADs = 1; skew-uniform slope)
#' @param intercept_speTADs intercept of a specific SKR used as a baseline (default: intercept_speTADs = 1.86; skew-uniform intercept)
#' @param lin_mod Indicates the type of linear model to use for (SKR): choose "lm" or "mblm"
#' @param distance_metric Indicates the method to compute distance-based regression parameters: choose "RMSE" (for Root Mean Square Error, default) or "MAE" (for Mean Absolute Error)
#' @returns RDS files with:
#' abundance for observed and randomized communities,
#' moments (mean, variance, skewness & kurtosis) for observed and randomized communities,
#' Standardized Effect Size (SES) values of the moments for observed compared to randomized communities and significance,
#' SKR parameters for observed and randomized communities
#' SES of the SKR parameters (slope, intercept, Rsquare, distance from predicted TADs, distance from specific TADs & CV of the distance from specific TADs) for observed compared to randomized communities and significance.
#' @export
#' @examples
#'
#' Example of function used, with "abundance" dataframe of grassland plant communities observed under different management practices over time, and "trait" dataframe functional trait per species  (SLA) from TRY database (Kattge et al. 2020)
#'
#' SKR.TAD::DataAnalysisTAD(
#' weights = SKR.TAD::abundance[,5:102],
#' weightsFactor = SKR.TAD::abundance[,c("Year", "Plot", "Treatment", "Bloc")],
#' dataToTreat = log(SKR.TAD::trait[["SLA"]]),
#' aggregationFactorName = c("Year", "Bloc"),
#' statisticsFactorName = c("Treatment"),
#' regenerateAbundanceDataFrame = T,
#' regenerateWeightedMomentsDataFrame = T,
#' regenerateStatPerObsDataFrame = T,
#' regenerateStatPerRandDataFrame = T,
#' randomizationNumber = 1000,
#' seed = 666,
#' abundanceDataFrameRDS = "./Output/abundanceDataFrame.RDS",
#' weightedMomentsDataFrameRDS = "./Output/MomentsDataFrame.RDS",
#' statPerObsDataFrameRDS = "./Output/SES_MomentsDataFrame.RDS",
#' statPerRandDataFrameRDS = "./Output/SKRDataFrame.RDS",
#' statSKRparam = "./Output/SES_SKRDataFrame.RDS",
#' significanceThreshold = c(0.05, 0.95),
#' slope_speTADs = 1,
#' intercept_speTADs = 1.86,
#' distance_metric = "RMSE",
#' lin_mod = "lm"
#' )

DataAnalysisTAD <- function(
    weights,
    weightsFactor,
    dataToTreat,
    randomizationNumber,
    aggregationFactorName = NULL,
    statisticsFactorName = NULL,
    seed = 123456,
    abundanceDataFrameRDS = NULL,
    weightedMomentsDataFrameRDS = NULL,
    statPerObsDataFrameRDS = NULL,
    statPerRandDataFrameRDS = NULL,
    statSKRparam,
    regenerateAbundanceDataFrame = FALSE,
    regenerateWeightedMomentsDataFrame = FALSE,
    regenerateStatPerObsDataFrame = FALSE,
    regenerateStatPerRandDataFrame = FALSE,
    significanceThreshold = c(0.05, 0.95),
    doParallel = TRUE,
    lin_mod = "lm",
    slope_speTADs = 1,
    intercept_speTADs = 1.86,
    distance_metric = "RMSE"
) {

  # preliminary test on input data
  if (nrow(weights) != nrow(weightsFactor)) {
    stop("weights and weightsFactor must have the same number of rows !")
  }

  if (ncol(weights) != length(dataToTreat)) {
    stop("the number of column of weights data must be equals to the length of the data to treat !")
  }

  # Generate or load random matrix
  if (is.null(abundanceDataFrameRDS) ||
      (!is.null(abundanceDataFrameRDS) && !file.exists(abundanceDataFrameRDS)) ||
      regenerateAbundanceDataFrame) {

    if (is.null(aggregationFactorName)) {
      aggregationFactor <- NULL
    }else {
      aggregationFactor <- as.data.frame(weightsFactor[, aggregationFactorName])
    }

    abundanceDataframe <- generateRandomMatrix(weights = weights,
                                               aggregationFactor = aggregationFactor,
                                               randomizationNumber = randomizationNumber,
                                               seed = seed,
                                               weightsDataframeRDS = abundanceDataFrameRDS)
  }else {
    abundanceDataframe <- readRDS(file = abundanceDataFrameRDS)
  }

  # Remove the species which have no trait value
  speciesToRemove <- which(is.na(dataToTreat))
  traitData <- dataToTreat

  if (length(speciesToRemove) != 0) {
    traitData <- traitData[-speciesToRemove]

    # Remove the species which have no trait value
    weights <- weights[, -speciesToRemove]
    abundanceDataframe <- abundanceDataframe[, -(1 + speciesToRemove)]
  }

  # Remove the observation with a total abundance of 0
  weightsToRemove <- which(rowSums(weights) == 0)

  if (length(weightsToRemove) != 0) {
    weightsFactor <- weightsFactor[-weightsToRemove, ]
    weights <- weights[-weightsToRemove, ]
  }

  abundanceToRemove <- which(rowSums(abundanceDataframe[, 2:ncol(abundanceDataframe)]) == 0)

  if (length(abundanceToRemove) != 0) {
    abundanceDataframe <- abundanceDataframe[-abundanceToRemove, ]
  }

  # Generate or load moments dataframe
  if (is.null(weightedMomentsDataFrameRDS) ||
      (!is.null(weightedMomentsDataFrameRDS) && !file.exists(weightedMomentsDataFrameRDS)) ||
      regenerateStatPerObsDataFrame) {
    # Compute for each line the weighted mean, variance, skewness, kurtosis and distance to lower boundary
    weightedMomentsList <-
      r4urep::weightedMVSK(data = traitData, weights = as.matrix(abundanceDataframe[, 2:(length(traitData) + 1)]))

    # Create a dataframe with the weighted moments and save it
    weightedMoments <- data.frame(matrix(data = NA, ncol = 0, nrow = nrow(abundanceDataframe)))
    weightedMoments$Number <- abundanceDataframe$Number
    weightedMoments$mean <- weightedMomentsList[["mean"]]
    weightedMoments$variance <- weightedMomentsList[["variance"]]
    weightedMoments$skewness <- weightedMomentsList[["skewness"]]
    weightedMoments$kurtosis <- weightedMomentsList[["kurtosis"]]
    rm(weightedMomentsList)
    weightedMoments$distance_speTADs <- weightedMoments$kurtosis - (slope_speTADs*weightedMoments$skewness*weightedMoments$skewness + intercept_speTADs)

    weightedMoments <-
      cbind(weightsFactor[rep(x = seq_len(nrow(weightsFactor)), times = randomizationNumber + 1), ,
                          drop = FALSE],
            weightedMoments)

    if (!is.null(weightedMomentsDataFrameRDS)) {
      saveRDS(object = weightedMoments, file = weightedMomentsDataFrameRDS)
    }
  } else {
    weightedMoments <- readRDS(file = weightedMomentsDataFrameRDS)
  }

  # Generate or load statistics per observation dataframe
  if (is.null(statPerObsDataFrameRDS) ||
      (!is.null(statPerObsDataFrameRDS) && !file.exists(statPerObsDataFrameRDS)) ||
      regenerateStatPerRandDataFrame) {
    # compute statistics for null model for mean/var/skewness/kurtosis
    statisticsPerObservation <- weightsFactor

    for (i in seq_len(nrow(statisticsPerObservation))){
      statisticsPerObservation[i, (ncol(weightsFactor) + 1):(ncol(weightsFactor) + 4)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = weightedMoments$mean[i],
          randomValues = weightedMoments$mean[(1:randomizationNumber) * nrow(statisticsPerObservation) + i],
          significanceThreshold = significanceThreshold)
      statisticsPerObservation[i, (ncol(weightsFactor) + 5):(ncol(weightsFactor) + 8)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = weightedMoments$variance[i],
          randomValues = weightedMoments$variance[(1:randomizationNumber) * nrow(statisticsPerObservation) + i],
          significanceThreshold = significanceThreshold)
      statisticsPerObservation[i, (ncol(weightsFactor) + 9):(ncol(weightsFactor) + 12)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = weightedMoments$skewness[i],
          randomValues = weightedMoments$skewness[(1:randomizationNumber) * nrow(statisticsPerObservation) + i],
          significanceThreshold = significanceThreshold)
      statisticsPerObservation[i, (ncol(weightsFactor) + 13):(ncol(weightsFactor) + 16)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = weightedMoments$kurtosis[i],
          randomValues = weightedMoments$kurtosis[(1:randomizationNumber) * nrow(statisticsPerObservation) + i],
          significanceThreshold = significanceThreshold)
    }
    commonColName <- c("standardizedObserved",
                       "standardizedMinQuantile",
                       "standardizedMaxQuantile",
                       "significance")
    colnames(statisticsPerObservation) <- c(colnames(weightsFactor),
                                            paste0(commonColName, "Mean"),
                                            paste0(commonColName, "Variance"),
                                            paste0(commonColName, "Skewness"),
                                            paste0(commonColName, "Kurtosis"))

    if (!is.null(statPerObsDataFrameRDS)) {
      saveRDS(object = statisticsPerObservation, file = statPerObsDataFrameRDS)
    }
  }

  # Generate or load statistics per observation dataframe
  if (is.null(statPerRandDataFrameRDS) ||
      (!is.null(statPerRandDataFrameRDS) && !file.exists(statPerRandDataFrameRDS)) ||
      regenerateWeightedMomentsDataFrame) {

    # Construct the id for statistics
    if (!is.null(statisticsFactorName)) {
      statisticsId <- apply(as.data.frame(weightsFactor[, statisticsFactorName]), 1, paste, collapse = "_")
    }else {
      statisticsId <- rep(x = "_", times = nrow(weightsFactor))
    }

    # Construct a list which contains for each statistics factor the species which are valid,
    # i.e. the sum of abundance is not equal to 0
    statisticsFactorSpeciesList <- list()
    for (statFactor in unique(statisticsId)) {
      statisticsFactorSpeciesList[[statFactor]] <-
        as.vector(which(colSums(weights[which(statisticsId == statFactor), ]) != 0))
    }

    # Generate the SKR analysis per null model regarding the factor given in parameter
    statisticsPerRandom <- data.frame(matrix(data = NA, nrow = (randomizationNumber + 1) *
                                               length(statisticsFactorSpeciesList), ncol = 0))
    lengthFactor <- length(names(statisticsFactorSpeciesList))
    abundanceDataframe$skewness <- weightedMoments$skewness
    abundanceDataframe$kurtosis <- weightedMoments$kurtosis
    abundanceDataframe$distance_speTADs <- weightedMoments$distance_speTADs

    for (i in 0:randomizationNumber) {
      for (j in 1:lengthFactor) {
        statisticsPerRandom$Number[i * lengthFactor + j] <- i

        statisticsPerRandom[i * lengthFactor + j, statisticsFactorName] <-
          weightsFactor[which(statisticsId == names(statisticsFactorSpeciesList)[j])[1], statisticsFactorName]

        dfToAnalyze <- abundanceDataframe[which(x = abundanceDataframe$Number == i), ]
        dfToAnalyze <- dfToAnalyze[which(x = statisticsId == names(statisticsFactorSpeciesList)[j]), ]
        y <- dfToAnalyze$kurtosis
        x <- dfToAnalyze$skewness^2
        distance_speTADs <- dfToAnalyze$distance_speTADs^2

        # for lintr
        if(lin_mod == "lm"){
          x <- x
          fit <- lm(y ~ x)
        }else if(lin_mod == "mblm"){
          x <- x
          fit <- mblm::mblm(y ~ x)
        }

        statisticsPerRandom$Slope[i * lengthFactor + j] <- fit$coefficients[2]
        statisticsPerRandom$Intercept[i * lengthFactor + j] <- fit$coefficients[1]
        statisticsPerRandom$Rsquare[i * lengthFactor + j] <-
          1 - (mean(stats::residuals(fit)^2, na.rm = TRUE) / stats::var(y, na.rm = TRUE))
        if(distance_metric == "RMSE"){
          statisticsPerRandom$distance_predicted_TADs[i * lengthFactor + j] <- sqrt(mean(fit$residuals^2, na.rm = TRUE))
          statisticsPerRandom$distance_specific_TADs[i * lengthFactor + j] <- sqrt(mean(distance_speTADs, na.rm = T))
        }else if (distance_metric == "MAE"){
          statisticsPerRandom$distance_predicted_TADs[i * lengthFactor + j] <- mean(sqrt(fit$residuals^2), na.rm = TRUE)
          statisticsPerRandom$distance_specific_TADs[i * lengthFactor + j] <- mean(sqrt(distance_speTADs), na.rm = T)
        }
        statisticsPerRandom$CV_distance_specific_TADs[i * lengthFactor + j] <- sd(sqrt(distance_speTADs), na.rm = T)*100/mean(sqrt(distance_speTADs), na.rm = T)
      }
    }

    if (!is.null(statPerRandDataFrameRDS)) {
      saveRDS(object = statisticsPerRandom, file = statPerRandDataFrameRDS)
    }
  }else {
    statisticsPerRandom <- readRDS(file = statPerRandDataFrameRDS)
  }
  # Compute SES of the SKR parameters
  SKRparam <- readRDS(file = statPerRandDataFrameRDS)
  SES_SKR <- data.frame()
  for (i in unique(SKRparam[[statisticsFactorName]])){
    SES_SKR <- rbind(SES_SKR,
                     data.frame(
                       Slope_SES = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Slope,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Slope,
                         significanceThreshold = significanceThreshold)[[1]][1],
                       Slope_Signi = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Slope,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Slope,
                         significanceThreshold = significanceThreshold)[[4]][1],
                       Intercept_SES = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Intercept,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Intercept,
                         significanceThreshold = significanceThreshold)[[1]][1],
                       Intercept_Signi = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Intercept,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Intercept,
                         significanceThreshold = significanceThreshold)[[4]][1],
                       Rsquare_SES = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Rsquare,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Rsquare,
                         significanceThreshold = significanceThreshold)[[1]][1],
                       Rsquare_Signi = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Rsquare,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Rsquare,
                         significanceThreshold = significanceThreshold)[[4]][1],
                       distance_predicted_TADs_SES = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$distance_predicted_TADs,
                         randomValues = SKRparam[SKRparam$Number > 0,]$distance_predicted_TADs,
                         significanceThreshold = significanceThreshold)[[1]][1],
                       distance_predicted_TADs_Signi = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$distance_predicted_TADs,
                         randomValues = SKRparam[SKRparam$Number > 0,]$distance_predicted_TADs,
                         significanceThreshold = significanceThreshold)[[4]][1],
                       distance_specific_TADs_SES = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$distance_specific_TADs,
                         randomValues = SKRparam[SKRparam$Number > 0,]$distance_specific_TADs,
                         significanceThreshold = significanceThreshold)[[1]][1],
                       distance_specific_TADs_Signi = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$distance_specific_TADs,
                         randomValues = SKRparam[SKRparam$Number > 0,]$distance_specific_TADs,
                         significanceThreshold = significanceThreshold)[[4]][1],
                       CV_distance_specific_TADs_SES = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$CV_distance_specific_TADs,
                         randomValues = SKRparam[SKRparam$Number > 0,]$CV_distance_specific_TADs,
                         significanceThreshold = significanceThreshold)[[1]][1],
                       CV_distance_specific_TADs_Signi = r4urep::nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$CV_distance_specific_TADs,
                         randomValues = SKRparam[SKRparam$Number > 0,]$CV_distance_specific_TADs,
                         significanceThreshold = significanceThreshold)[[4]][1],
                       statisticsFactorName = i))
  }
  names(SES_SKR)[names(SES_SKR) == "statisticsFactorName"] <- statisticsFactorName

  saveRDS(object = SES_SKR, file = statSKRparam)
}
