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
#' @param path_MomentsDataFrame the path and name of the RDS file to load/save the dataframe which contains the calculated moments (for observations and randomizations)
#' @param path_SESMomentsDataFrame the path and name of the RDS file to load/save the dataframe which contains the moments standardized effect size (observations regarding randomizations)
#' @param path_SKRDataFrame the path and name of the RDS file to load/save the dataframe which contains the SKR parameters (for observations and randomizations)
#' @param path_SESSKRDataFrame the path and name of the RDS file to load/save the dataframe which contains the SKR parameters standardized effect size (observations regarding randomizations)
#' @param regenerate_abundanceDataFrame boolean to specify if the abundance dataframe is computed again
#' @param regenerate_MomentsDataFrame boolean to specify if the moments dataframe is computed again
#' @param regenerate_SESMomentsDataFrame boolean to specify if the moments standardized effect size dataframe is computed again
#' @param regenerate_SKRDataFrame boolean to specify if the SKR parameters dataframe is computed again
#' @param regenerate_SESSKRDataFrame boolean to specify if the SKR parameters standardized effect size dataframe is computed again
#' @param significanceThreshold the significance threshold to consider that the observed value is in the randomized value
#' @param doParallel Indicates if we use parallelism to construct the random matrix
#' @param slope_ref_TADs slope of a reference SKR used as a baseline (default: slope_ref_TADs = 1; skew-uniform slope)
#' @param intercept_ref_TADs intercept of a reference SKR used as a baseline (default: intercept_ref_TADs = 1.86; skew-uniform intercept)
#' @param lin_mod Indicates the type of linear model to use for (SKR): choose "lm" or "mblm"
#' @param distance_metric Indicates the method to compute distance-based regression parameters: choose "RMSE" (for Root Mean Square Error, default) or "MAE" (for Mean Absolute Error)
#' @param log_distance Logical. If `TRUE`, the error metrics (RMSE, MAE) are computed on the log-transformed distances between observed and predicted kurtosis values.
#' @returns RDS files with:
#' abundance for observed and randomized communities,
#' moments (mean, variance, skewness & kurtosis) for observed and randomized communities,
#' Standardized Effect Size (SES) values of the moments for observed compared to randomized communities and significance,
#' SKR parameters for observed and randomized communities
#' SES of the SKR parameters (slope, intercept, Rsquare, distance from predicted TADs, distance from reference TADs & CV of the distance from reference TADs).
#' @export
#' @examples
#' head(abundance)
#' head(trait)
#' DataAnalysisTAD(
#' Abundance = abundance[,5:102],
#' AbundanceFactor = abundance[,c("Year", "Plot", "Treatment", "Bloc")],
#' TraitData = log(trait[["SLA"]] + 1),
#' randomizationFactorName = c("Year", "Bloc"),
#' statisticsFactorName = c("Treatment"),
#' regenerate_abundanceDataFrame = T,
#' regenerate_MomentsDataFrame = T,
#' regenerate_SESMomentsDataFrame = T,
#' regenerate_SKRDataFrame = T,
#' regenerate_SESSKRDataFrame = T,
#' randomizationNumber = 10,
#' seed = 666,
#' path_abundanceDataFrame = NULL,
#' path_MomentsDataFrame = NULL,
#' path_SESMomentsDataFrame = NULL,
#' path_SKRDataFrame = NULL,
#' path_SESSKRDataFrame = NULL,
#' significanceThreshold = c(0.05, 0.95),
#' slope_ref_TADs = 1,
#' intercept_ref_TADs = 1.86,
#' distance_metric = "RMSE",
#' lin_mod = "lm",
#' log_distance = FALSE,
#' doParallel = FALSE
#' )

DataAnalysisTAD <- function(
    Abundance,
    AbundanceFactor,
    TraitData,
    randomizationNumber,
    randomizationFactorName = NULL,
    statisticsFactorName = NULL,
    seed = 123456,
    path_abundanceDataFrame = NULL,
    path_MomentsDataFrame = NULL,
    path_SESMomentsDataFrame = NULL,
    path_SKRDataFrame = NULL,
    path_SESSKRDataFrame = NULL,
    regenerate_abundanceDataFrame = FALSE,
    regenerate_MomentsDataFrame = FALSE,
    regenerate_SESMomentsDataFrame = FALSE,
    regenerate_SKRDataFrame = FALSE,
    regenerate_SESSKRDataFrame = FALSE,
    significanceThreshold = c(0.05, 0.95),
    doParallel = TRUE,
    lin_mod = "lm",
    slope_ref_TADs = 1,
    intercept_ref_TADs = 1.86,
    distance_metric = "RMSE",
    log_distance = FALSE
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

    abundanceDataframe <- AbundanceRandomization(Abundance = Abundance,
                                                 randomizationFactor = randomizationFactor,
                                                 randomizationNumber = randomizationNumber,
                                                 seed = seed,
                                                 path_abundanceDataFrame = path_abundanceDataFrame,
                                                 doParallel = doParallel)
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
  if (is.null(path_MomentsDataFrame) ||
      (!is.null(path_MomentsDataFrame) && !file.exists(path_MomentsDataFrame)) ||
      regenerate_MomentsDataFrame) {
    # Compute for each line the weighted mean, variance, skewness, kurtosis and distance to reference TADs (with the mention slope_ref_TADs and intercept_ref_TADs)
    MomentsList <-
      r4urep::weightedMVSK(data = traitData, weights = as.matrix(abundanceDataframe[, 2:(length(traitData) + 1)]))

    # Create a dataframe with the weighted moments and save it
    MomentsDataFrame <- data.frame(matrix(data = NA, ncol = 0, nrow = nrow(abundanceDataframe)))
    MomentsDataFrame$Number <- abundanceDataframe$Number
    MomentsDataFrame$mean <- MomentsList[["mean"]]
    MomentsDataFrame$variance <- MomentsList[["variance"]]
    MomentsDataFrame$skewness <- MomentsList[["skewness"]]
    MomentsDataFrame$kurtosis <- MomentsList[["kurtosis"]]
    rm(MomentsList)
    MomentsDataFrame$distance_reference_TADs <- MomentsDataFrame$kurtosis - (slope_ref_TADs*MomentsDataFrame$skewness^2 + intercept_ref_TADs)

    MomentsDataFrame <-
      cbind(AbundanceFactor[rep(x = seq_len(nrow(AbundanceFactor)), times = randomizationNumber + 1), ,
                            drop = FALSE],
            MomentsDataFrame)

    if (!is.null(path_MomentsDataFrame)) {
      message("Save MomentsDataFrame", path_MomentsDataFrame)
      saveRDS(object = MomentsDataFrame, file = path_MomentsDataFrame)
    }
  } else {
    MomentsDataFrame <- readRDS(file = path_MomentsDataFrame)
  }

  # Generate or load moments standardized effect size (SES)
  if (is.null(path_SESMomentsDataFrame) ||
      (!is.null(path_SESMomentsDataFrame) && !file.exists(path_SESMomentsDataFrame)) ||
      regenerate_SESMomentsDataFrame) {

    # compute moments standardized effect size (SES of mean/var/skewness/kurtosis)
    SESMomentsDataFrame <- AbundanceFactor

    for (i in seq_len(nrow(SESMomentsDataFrame))){
      SESMomentsDataFrame[i, (ncol(AbundanceFactor) + 1):(ncol(AbundanceFactor) + 4)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = MomentsDataFrame$mean[i],
          randomValues = MomentsDataFrame$mean[(1:randomizationNumber) * nrow(SESMomentsDataFrame) + i],
          significanceThreshold = significanceThreshold)
      SESMomentsDataFrame[i, (ncol(AbundanceFactor) + 5):(ncol(AbundanceFactor) + 8)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = MomentsDataFrame$variance[i],
          randomValues = MomentsDataFrame$variance[(1:randomizationNumber) * nrow(SESMomentsDataFrame) + i],
          significanceThreshold = significanceThreshold)
      SESMomentsDataFrame[i, (ncol(AbundanceFactor) + 9):(ncol(AbundanceFactor) + 12)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = MomentsDataFrame$skewness[i],
          randomValues = MomentsDataFrame$skewness[(1:randomizationNumber) * nrow(SESMomentsDataFrame) + i],
          significanceThreshold = significanceThreshold)
      SESMomentsDataFrame[i, (ncol(AbundanceFactor) + 13):(ncol(AbundanceFactor) + 16)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = MomentsDataFrame$kurtosis[i],
          randomValues = MomentsDataFrame$kurtosis[(1:randomizationNumber) * nrow(SESMomentsDataFrame) + i],
          significanceThreshold = significanceThreshold)
    }

    commonColName <- c("SES",
                       "SESMinQuantile",
                       "SESMaxQuantile",
                       "significance")

    colnames(SESMomentsDataFrame) <- c(colnames(AbundanceFactor),
                                       paste0(commonColName, "mean"),
                                       paste0(commonColName, "variance"),
                                       paste0(commonColName, "skewness"),
                                       paste0(commonColName, "kurtosis"))

    if (!is.null(path_SESMomentsDataFrame)) {
      message("Save SESMomentsDataFrame", path_SESMomentsDataFrame)
      saveRDS(object = SESMomentsDataFrame, file = path_SESMomentsDataFrame)
    }
  } else {
    SESMomentsDataFrame <- readRDS(file = path_SESMomentsDataFrame)
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
    SKRDataFrame <- data.frame(matrix(data = NA, nrow = (randomizationNumber + 1) *
                                        length(statisticsFactorSpeciesList), ncol = 0))
    lengthFactor <- length(names(statisticsFactorSpeciesList))
    abundanceDataframe$skewness <- MomentsDataFrame$skewness
    abundanceDataframe$kurtosis <- MomentsDataFrame$kurtosis
    abundanceDataframe$distance_reference_TADs <- MomentsDataFrame$distance_reference_TADs

    for (i in 0:randomizationNumber) {
      for (j in 1:lengthFactor) {
        SKRDataFrame$Number[i * lengthFactor + j] <- i

        SKRDataFrame[i * lengthFactor + j, statisticsFactorName] <-
          AbundanceFactor[which(statisticsId == names(statisticsFactorSpeciesList)[j])[1], statisticsFactorName]

        dfToAnalyze <- abundanceDataframe[which(x = abundanceDataframe$Number == i), ]
        dfToAnalyze <- dfToAnalyze[which(x = statisticsId == names(statisticsFactorSpeciesList)[j]), ]
        y <- dfToAnalyze$kurtosis
        x <- dfToAnalyze$skewness^2

        # for lintr
        if(lin_mod == "lm"){
          x <- x
          fit <- stats::lm(y ~ x)
        }else if(lin_mod == "mblm"){
          x <- x
          fit <- mblm::mblm(y ~ x)
        }

        y_pred <- fitted(fit)
        y_pred_ref <- slope_ref_TADs * x + intercept_ref_TADs
        residuals_fit <- stats::residuals(fit)

        SKRDataFrame$Slope[i * lengthFactor + j] <- fit$coefficients[2]
        SKRDataFrame$Intercept[i * lengthFactor + j] <- fit$coefficients[1]
        SKRDataFrame$Rsquare[i * lengthFactor + j] <-
          1 - (mean(residuals_fit^2, na.rm = TRUE) / stats::var(y, na.rm = TRUE))

        if (distance_metric == "RMSE") {
          if (!log_distance) {
            SKRDataFrame$distance_predicted_TADs[i * lengthFactor + j] <- sqrt(mean(residuals_fit^2, na.rm = TRUE))
            SKRDataFrame$distance_reference_TADs[i * lengthFactor + j] <- sqrt(mean((y - y_pred_ref)^2, na.rm = TRUE))
          } else {
            SKRDataFrame$distance_predicted_TADs[i * lengthFactor + j] <- sqrt(mean((log(y + 1) - log(y_pred + 1))^2, na.rm = TRUE))
            SKRDataFrame$distance_reference_TADs[i * lengthFactor + j] <- sqrt(mean((log(y + 1) - log(y_pred_ref + 1))^2, na.rm = TRUE))
          }
        } else if (distance_metric == "MAE") {
          if (!log_distance) {
            SKRDataFrame$distance_predicted_TADs[i * lengthFactor + j] <- mean(abs(residuals_fit), na.rm = TRUE)
            SKRDataFrame$distance_reference_TADs[i * lengthFactor + j] <- mean(abs(y - y_pred_ref), na.rm = TRUE)
          } else {
            SKRDataFrame$distance_predicted_TADs[i * lengthFactor + j] <- mean(abs(log(y + 1) - log(y_pred + 1)), na.rm = TRUE)
            SKRDataFrame$distance_reference_TADs[i * lengthFactor + j] <- mean(abs(log(y + 1) - log(y_pred_ref + 1)), na.rm = TRUE)
          }
        }

        SKRDataFrame$CV_distance_reference_TADs[i * lengthFactor + j] <- sd(abs(y - y_pred_ref), na.rm = TRUE) * 100 / mean(abs(y - y_pred_ref), na.rm = TRUE)
      }
    }

    if (!is.null(path_SKRDataFrame)) {
      message("Save SKRDataFrame", path_SKRDataFrame)
      saveRDS(object = SKRDataFrame, file = path_SKRDataFrame)
    }
  } else {
    SKRDataFrame <- readRDS(file = path_SKRDataFrame)
  }

  # Generate or load statistics on SKR
  if (is.null(path_SESSKRDataFrame) ||
      (!is.null(path_SESSKRDataFrame) && !file.exists(path_SESSKRDataFrame)) ||
      regenerate_SESSKRDataFrame) {

    # compute statistics (SES) for null model for Slope, Intercept, Rsquare, distance_predicted_TADs, distance_reference_TADs, CV_distance_reference_TADs
    SESSKRDataFrame <- unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))

    for (i in seq_len(nrow(SESSKRDataFrame))){
      SESSKRDataFrame[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 1):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 4)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRDataFrame$Slope[i],
          randomValues = SKRDataFrame$Slope[(1:randomizationNumber) * nrow(SESSKRDataFrame) + i],
          significanceThreshold = significanceThreshold)
      SESSKRDataFrame[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 5):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 8)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRDataFrame$Intercept[i],
          randomValues = SKRDataFrame$Intercept[(1:randomizationNumber) * nrow(SESSKRDataFrame) + i],
          significanceThreshold = significanceThreshold)
      SESSKRDataFrame[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 9):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 12)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRDataFrame$Rsquare[i],
          randomValues = SKRDataFrame$Rsquare[(1:randomizationNumber) * nrow(SESSKRDataFrame) + i],
          significanceThreshold = significanceThreshold)
      SESSKRDataFrame[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 13):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 16)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRDataFrame$distance_predicted_TADs[i],
          randomValues = SKRDataFrame$distance_predicted_TADs[(1:randomizationNumber) * nrow(SESSKRDataFrame) + i],
          significanceThreshold = significanceThreshold)
      SESSKRDataFrame[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 17):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 20)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRDataFrame$distance_reference_TADs[i],
          randomValues = SKRDataFrame$distance_reference_TADs[(1:randomizationNumber) * nrow(SESSKRDataFrame) + i],
          significanceThreshold = significanceThreshold)
      SESSKRDataFrame[i, (ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 21):(ncol(unique(as.data.frame(AbundanceFactor[, statisticsFactorName]))) + 24)] <-
        r4urep::nullModelDistributionStatistics(
          observedValue = SKRDataFrame$CV_distance_reference_TADs[i],
          randomValues = SKRDataFrame$CV_distance_reference_TADs[(1:randomizationNumber) * nrow(SESSKRDataFrame) + i],
          significanceThreshold = significanceThreshold)
    }
    commonColName <- c("SES",
                       "SESMinQuantile",
                       "SESMaxQuantile",
                       "significance")

    colnames(SESSKRDataFrame) <- c(paste0(statisticsFactorName),
                                   paste0(commonColName, "Slope"),
                                   paste0(commonColName, "Intercept"),
                                   paste0(commonColName, "Rsquare"),
                                   paste0(commonColName, "distance_predicted_TADs"),
                                   paste0(commonColName, "distance_reference_TADs"),
                                   paste0(commonColName, "CV_distance_reference_TADs"))

    if (!is.null(path_SESSKRDataFrame)) {
      message("Save SESSKRDataFrame", path_SESSKRDataFrame)
      saveRDS(object = SESSKRDataFrame, file = path_SESSKRDataFrame)
    }
  } else {
    SESSKRDataFrame <- readRDS(file = path_SESSKRDataFrame)
  }

  return(list(
    abundanceDataFrame = abundanceDataFrame,
    MomentsDataFrame = MomentsDataFrame,
    SESMomentsDataFrame = SESMomentsDataFrame,
    SKRDataFrame = SKRDataFrame,
    SESSKRDataFrame = SESSKRDataFrame
  ))

  list(
    abundanceDataFrame = abundanceDataFrame,
    MomentsDataFrame = MomentsDataFrame,
    SESMomentsDataFrame = SESMomentsDataFrame,
    SKRDataFrame = SKRDataFrame,
    SESSKRDataFrame = SESSKRDataFrame
  )
}
