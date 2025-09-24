library(devtools) #devtools_2.4.5
devtools::install_gitlab(repo = "urep/dev_utils/r_utils/r4urep",
                         host = "https://forge.inrae.fr",
                         force = T)

## Step 2. Function: SKR ANALYSIS ----

#' @title Launch the analysis of the TADs with the SKR framework
#' @description Launch the SKR analysis of the TADs, and generate output dataset (distribution moments, SKR parameters, and statistics)
#' @param Abundance the dataframe of abundance (or related weights measure), one row correspond to a series of observation
#' @param AbundanceFactor the dataframe which contains the different factor linked to the abundance dataframe
#' @param TraitData a vector of the trait data linked to the different abundance factor
#' @param randomizationNumber the number of random abundance matrix to generate
#' @param randomizationFactorName vector of factor name for the generation of random matrix
#' @param statisticsFactorName vector of factor name for the computation of statistics for each generated matrix
#' @param seed the seed of the pseudo random number generator
#' @param path_abundanceDataFrame the path and name of the RDS file to load/save the dataframe which contains the observed abundance data and the generated random abundance matrix
#' @param path_momentsDataFrame the path and name of the RDS file to load/save the dataframe which contains the calculated moments (for observations and randomizations)
#' @param path_SESmomentsDataFrame the path and name of the RDS file to load/save the dataframe which contains the moments standardized effect size (observations regarding randomizations)
#' @param path_SKRDataFrame the path and name of the RDS file to load/save the dataframe which contains the SKR parameters (for observations and randomizations)
#' @param path_SESSKRDataFrame the path and name of the RDS file to load/save the dataframe which contains the SKR parameters standardized effect size (observations regarding randomizations)
#' @param regenerate_abundanceDataFrame boolean to specify if the abundance dataframe is computed again
#' @param regenerate_momentsDataFrame boolean to specify if the moments dataframe is computed again
#' @param regenerate_SESmomentsDataFrame boolean to specify if the moments standardized effect size dataframe is computed again
#' @param regenerate_SKRDataFrame boolean to specify if the SKR parameters dataframe is computed again
#' @param regenerate_SESSKRDataFrame boolean to specify if the SKR parameters standardized effect size dataframe is computed again
#' @param significanceThreshold the significance threshold to consider that the observed value is in the randomized value
#' @param doParallel Indicates if we use parallelism to construct the random matrix
#' @param slope_ref_TADs slope of a reference SKR used as a baseline (default: slope_ref_TADs = 1; skew-uniform slope)
#' @param intercept_ref_TADs intercept of a reference SKR used as a baseline (default: intercept_ref_TADs = 1.86; skew-uniform intercept)
#' @param lin_mod Indicates the type of linear model to use for (SKR): choose "lm" or "mblm"
#' @param distance_metric Indicates the method to compute distance-based regression parameters: choose "RMSE" (for Root Mean Square Error, default) or "MAE" (for Mean Absolute Error)
#' @returns RDS files with:
#' abundance for observed and randomized communities,
#' moments (mean, variance, skewness & kurtosis) for observed and randomized communities,
#' Standardized Effect Size (SES) values of the moments for observed compared to randomized communities and significance,
#' SKR parameters for observed and randomized communities
#' SES of the SKR parameters (slope, intercept, Rsquare, distance from predicted TADs, distance from reference TADs & CV of the distance from reference TADs).
#' @export
#' @examples
#' Example of how to use the function, with abundance dataframe of grassland
#' plant communities observed under contrasting management practices and
#' trait dataframe functional trait per species  (SLA) from TRY database (Kattge et al. 2020)
#' DataAnalysisTAD(Abundance = SKR.TAD::abundance[,5:102],
#' AbundanceFactor = SKR.TAD::abundance[,c("Year", "Plot", "Treatment", "Bloc")],
#' TraitData = log(SKR.TAD::trait[["SLA"]] + 1),
#' randomizationFactorName = c("Year", "Bloc"),
#' statisticsFactorName = c("Treatment"),
#' regenerate_abundanceDataFrame = T,
#' regenerate_momentsDataFrame = T,
#' regenerate_SESmomentsDataFrame = T,
#' regenerate_SKRDataFrame = T,
#' regenerate_SESSKRDataFrame = T,
#' randomizationNumber = 100,
#' seed = 666,
#' path_abundanceDataFrame = "./abundanceDataFrame.RDS",
#' path_momentsDataFrame = "./MomentsDataFrame.RDS",
#' path_SESmomentsDataFrame = "./SES_MomentsDataFrame.RDS",
#' path_SKRDataFrame = "./SKRDataFrame.RDS",
#' path_SESSKRDataFrame = "./SES_SKRDataFrame.RDS",
#' significanceThreshold = c(0.05, 0.95),
#' slope_ref_TADs = 1,
#' intercept_ref_TADs = 1.86,
#' distance_metric = "RMSE",
#' lin_mod = "lm")

DataAnalysisTAD <- function(
    Abundance,
    AbundanceFactor,
    TraitData,
    randomizationNumber,
    randomizationFactorName = NULL,
    statisticsFactorName = NULL,
    seed = 123456,
    path_abundanceDataFrame = NULL,
    path_momentsDataFrame = NULL,
    path_SESmomentsDataFrame = NULL,
    path_SKRDataFrame = NULL,
    path_SESSKRDataFrame,
    regenerate_abundanceDataFrame = FALSE,
    regenerate_momentsDataFrame = FALSE,
    regenerate_SESmomentsDataFrame = FALSE,
    regenerate_SKRDataFrame = FALSE,
    regenerate_SESSKRDataFrame = FALSE,
    significanceThreshold = c(0.05, 0.95),
    doParallel = TRUE,
    lin_mod = "lm",
    slope_ref_TADs = 1,
    intercept_ref_TADs = 1.86,
    distance_metric = "RMSE"
) {

  # preliminary test on input data
  if (nrow(Abundance) != nrow(AbundanceFactor)) {
    stop("Abundance and AbundanceFactor must have the same number of rows !")
  }

  if (ncol(Abundance) != length(TraitData)) {
    stop("the number of column of Abundance data must be equals to the length of the Trait data !")
  }

  # Generate or load random abundance matrix
  if (is.null(path_abundanceDataFrame) ||
      (!is.null(path_abundanceDataFrame) && !file.exists(path_abundanceDataFrame)) ||
      regenerate_abundanceDataFrame) {

    if (is.null(randomizationFactorName)) {
      randomizationFactor <- NULL
    }else {
      randomizationFactor <- as.data.frame(AbundanceFactor[, randomizationFactorName])
    }

    abundanceDataframe <- generateRandomMatrix(Abundance = Abundance,
                                               randomizationFactor = randomizationFactor,
                                               randomizationNumber = randomizationNumber,
                                               seed = seed,
                                               path_abundanceDataFrame = path_abundanceDataFrame)
  }else {
    abundanceDataframe <- readRDS(file = path_abundanceDataFrame)
  }

  # Remove the species which have no trait value
  speciesToRemove <- which(is.na(TraitData))
  traitData <- TraitData

  if (length(speciesToRemove) != 0) {
    traitData <- traitData[-speciesToRemove]

    # Remove the species which have no trait value
    Abundance <- Abundance[, -speciesToRemove]
    abundanceDataframe <- abundanceDataframe[, -(1 + speciesToRemove)]
  }

  # Remove the observation with a total abundance of 0
  AbundanceToRemove <- which(rowSums(Abundance) == 0)

  if (length(AbundanceToRemove) != 0) {
    AbundanceFactor <- AbundanceFactor[-AbundanceToRemove, ]
    Abundance <- Abundance[-AbundanceToRemove, ]
  }

  abundanceDataframeToRemove <- which(rowSums(abundanceDataframe[, 2:ncol(abundanceDataframe)]) == 0)

  if (length(abundanceDataframeToRemove) != 0) {
    abundanceDataframe <- abundanceDataframe[-abundanceDataframeToRemove, ]
  }

  # Generate or load moments dataframe
  if (is.null(path_momentsDataFrame) ||
      (!is.null(path_momentsDataFrame) && !file.exists(path_momentsDataFrame)) ||
      regenerate_momentsDataFrame) {
    # Compute for each line the weighted mean, variance, skewness, kurtosis and distance to reference TADs (with the mention slope_ref_TADs and intercept_ref_TADs)
    MomentsList <-
      r4urep::weightedMVSK(data = traitData, weights = as.matrix(abundanceDataframe[, 2:(length(traitData) + 1)]))

    # Create a dataframe with the weighted moments and save it
    Moments <- data.frame(matrix(data = NA, ncol = 0, nrow = nrow(abundanceDataframe)))
    Moments$Number <- abundanceDataframe$Number
    Moments$mean <- MomentsList[["mean"]]
    Moments$variance <- MomentsList[["variance"]]
    Moments$skewness <- MomentsList[["skewness"]]
    Moments$kurtosis <- MomentsList[["kurtosis"]]
    rm(MomentsList)
    Moments$distance_reference_TADs <- Moments$kurtosis - (slope_ref_TADs*Moments$skewness*Moments$skewness + intercept_ref_TADs)

    Moments <-
      cbind(AbundanceFactor[rep(x = seq_len(nrow(AbundanceFactor)), times = randomizationNumber + 1), ,
                            drop = FALSE],
            Moments)

    if (!is.null(path_momentsDataFrame)) {
      saveRDS(object = Moments, file = path_momentsDataFrame)
    }
  } else {
    Moments <- readRDS(file = path_momentsDataFrame)
  }

  # Generate or load moments standardized effect size (SES)
  if (is.null(path_SESmomentsDataFrame) ||
      (!is.null(path_SESmomentsDataFrame) && !file.exists(path_SESmomentsDataFrame)) ||
      regenerate_SESmomentsDataFrame) {

      # compute moments standardized effect size (SES of mean/var/skewness/kurtosis)
    SESMoments <- AbundanceFactor

    for (i in seq_len(nrow(SESMoments))){
      SESMoments[i, (ncol(AbundanceFactor) + 1):(ncol(AbundanceFactor) + 4)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = Moments$mean[i],
          randomValues = Moments$mean[(1:randomizationNumber) * nrow(SESMoments) + i],
          significanceThreshold = significanceThreshold)
      SESMoments[i, (ncol(AbundanceFactor) + 5):(ncol(AbundanceFactor) + 8)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = Moments$variance[i],
          randomValues = Moments$variance[(1:randomizationNumber) * nrow(SESMoments) + i],
          significanceThreshold = significanceThreshold)
      SESMoments[i, (ncol(AbundanceFactor) + 9):(ncol(AbundanceFactor) + 12)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = Moments$skewness[i],
          randomValues = Moments$skewness[(1:randomizationNumber) * nrow(SESMoments) + i],
          significanceThreshold = significanceThreshold)
      SESMoments[i, (ncol(AbundanceFactor) + 13):(ncol(AbundanceFactor) + 16)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = Moments$kurtosis[i],
          randomValues = Moments$kurtosis[(1:randomizationNumber) * nrow(SESMoments) + i],
          significanceThreshold = significanceThreshold)
    }

    commonColName <- c("SES",
                       "SESMinQuantile",
                       "SESMaxQuantile",
                       "significance")

    colnames(SESMoments) <- c(colnames(AbundanceFactor),
                                            paste0(commonColName, "Mean"),
                                            paste0(commonColName, "Variance"),
                                            paste0(commonColName, "Skewness"),
                                            paste0(commonColName, "Kurtosis"))

    if (!is.null(path_SESmomentsDataFrame)) {
      saveRDS(object = SESMoments, file = path_SESmomentsDataFrame)
    }
  }

  # Generate or load SKR parameters dataframe
  if (is.null(path_SKRDataFrame) ||
      (!is.null(path_SKRDataFrame) && !file.exists(path_SKRDataFrame)) ||
      regenerate_SKRDataFrame) {

    # Construct the id for statistics
    if (!is.null(statisticsFactorName)) {
      statisticsId <- apply(as.data.frame(AbundanceFactor[, statisticsFactorName]), 1, paste, collapse = "_")
    }else {
      statisticsId <- rep(x = "_", times = nrow(AbundanceFactor))
    }

    # Construct a list which contains for each statistics factor the species which are valid,
    # i.e. the sum of abundance is not equal to 0
    statisticsFactorSpeciesList <- list()
    for (statFactor in unique(statisticsId)) {
      statisticsFactorSpeciesList[[statFactor]] <-
        as.vector(which(colSums(Abundance[which(statisticsId == statFactor), ]) != 0))
    }

    # Generate the SKR analysis per null model regarding the factor given in parameter
    SKRparam <- data.frame(matrix(data = NA, nrow = (randomizationNumber + 1) *
                                               length(statisticsFactorSpeciesList), ncol = 0))
    lengthFactor <- length(names(statisticsFactorSpeciesList))
    abundanceDataframe$skewness <- Moments$skewness
    abundanceDataframe$kurtosis <- Moments$kurtosis
    abundanceDataframe$distance_reference_TADs <- Moments$distance_reference_TADs

    for (i in 0:randomizationNumber) {
      for (j in 1:lengthFactor) {
        SKRparam$Number[i * lengthFactor + j] <- i

        SKRparam[i * lengthFactor + j, statisticsFactorName] <-
          AbundanceFactor[which(statisticsId == names(statisticsFactorSpeciesList)[j])[1], statisticsFactorName]

        dfToAnalyze <- abundanceDataframe[which(x = abundanceDataframe$Number == i), ]
        dfToAnalyze <- dfToAnalyze[which(x = statisticsId == names(statisticsFactorSpeciesList)[j]), ]
        y <- dfToAnalyze$kurtosis
        x <- dfToAnalyze$skewness^2
        distance_reference_TADs <- dfToAnalyze$distance_reference_TADs^2

        # for lintr
        if(lin_mod == "lm"){
          x <- x
          fit <- stats::lm(y ~ x)
        }else if(lin_mod == "mblm"){
          x <- x
          fit <- mblm::mblm(y ~ x)
        }

        SKRparam$Slope[i * lengthFactor + j] <- fit$coefficients[2]
        SKRparam$Intercept[i * lengthFactor + j] <- fit$coefficients[1]
        SKRparam$Rsquare[i * lengthFactor + j] <-
          1 - (mean(stats::residuals(fit)^2, na.rm = TRUE) / stats::var(y, na.rm = TRUE))
        if(distance_metric == "RMSE"){
          SKRparam$distance_predicted_TADs[i * lengthFactor + j] <- sqrt(mean(fit$residuals^2, na.rm = TRUE))
          SKRparam$distance_reference_TADs[i * lengthFactor + j] <- sqrt(mean(distance_reference_TADs, na.rm = T))
        }else if (distance_metric == "MAE"){
          SKRparam$distance_predicted_TADs[i * lengthFactor + j] <- mean(sqrt(fit$residuals^2), na.rm = TRUE)
          SKRparam$distance_reference_TADs[i * lengthFactor + j] <- mean(sqrt(distance_reference_TADs), na.rm = T)
        }
        SKRparam$CV_distance_reference_TADs[i * lengthFactor + j] <- stats::sd(sqrt(distance_reference_TADs), na.rm = T)*100/mean(sqrt(distance_reference_TADs), na.rm = T)
      }
    }

    if (!is.null(path_SKRDataFrame)) {
      saveRDS(object = SKRparam, file = path_SKRDataFrame)
    }
  }else {
    SKRparam <- readRDS(file = path_SKRDataFrame)
  }

  # Generate or load statistics on SKR
  if (is.null(path_SESSKRDataFrame) ||
      (!is.null(path_SESSKRDataFrame) && !file.exists(path_SESSKRDataFrame)) ||
      regenerate_SESSKRDataFrame) {

    # compute statistics (SES) for null model for Slope, Intercept, Rsquare, distance_predicted_TADs, distance_reference_TADs, CV_distance_reference_TADs
    SES_SKRparam <- unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))

    for (i in seq_len(nrow(SES_SKRparam))){
      SES_SKRparam[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 1):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 4)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRparam$Slope[i],
          randomValues = SKRparam$Slope[(1:randomizationNumber) * nrow(SES_SKRparam) + i],
          significanceThreshold = significanceThreshold)
      SES_SKRparam[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 5):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 8)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRparam$Intercept[i],
          randomValues = SKRparam$Intercept[(1:randomizationNumber) * nrow(SES_SKRparam) + i],
          significanceThreshold = significanceThreshold)
      SES_SKRparam[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 9):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 12)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRparam$Rsquare[i],
          randomValues = SKRparam$Rsquare[(1:randomizationNumber) * nrow(SES_SKRparam) + i],
          significanceThreshold = significanceThreshold)
      SES_SKRparam[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 13):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 16)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRparam$distance_predicted_TADs[i],
          randomValues = SKRparam$distance_predicted_TADs[(1:randomizationNumber) * nrow(SES_SKRparam) + i],
          significanceThreshold = significanceThreshold)
      SES_SKRparam[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 17):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 20)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRparam$distance_reference_TADs[i],
          randomValues = SKRparam$distance_reference_TADs[(1:randomizationNumber) * nrow(SES_SKRparam) + i],
          significanceThreshold = significanceThreshold)
      SES_SKRparam[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 21):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 24)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRparam$CV_distance_reference_TADs[i],
          randomValues = SKRparam$CV_distance_reference_TADs[(1:randomizationNumber) * nrow(SES_SKRparam) + i],
          significanceThreshold = significanceThreshold)
    }
    commonColName <- c("SES",
                       "SESMinQuantile",
                       "SESMaxQuantile",
                       "significance")

    colnames(SES_SKRparam) <- c(paste0(statisticsFactorName),
                           paste0(commonColName, "Slope"),
                           paste0(commonColName, "Intercept"),
                           paste0(commonColName, "Rsquare"),
                           paste0(commonColName, "distance_predicted_TADs"),
                           paste0(commonColName, "distance_reference_TADs"),
                           paste0(commonColName, "CV_distance_reference_TADs"))

    if (!is.null(path_SESSKRDataFrame)) {
      saveRDS(object = SES_SKRparam, file = path_SESSKRDataFrame)
    }
  }
}
