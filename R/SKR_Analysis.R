# I. PACKAGES ----

library(devtools) #devtools_2.4.5
devtools::install_gitlab(repo = "urep/dev_utils/r_utils/r4urep",
                         host = "https://forge.inrae.fr")
library(parallel) #parallel_4.3.1
library(doParallel) #doParallel_1.0.17
library(foreach) #foreach_1.5.2
library(mblm) #mblm_0.12.1 
library(dplyr) #dplyr_1.1.3
library(ggplot2) #ggplot2_3.5.1
library(ggpubr) #ggpubr_0.6.0

# II. Function ----
## Step 1. FUNCTION: RANDOMIZATION ----
### a. Generate random matrix ----

#' @title Generate random matrix
#' @description Generate and save random matrix
#' @concept tad
#' @param weights the dataframe of weights, one row correspond to a series of observation
#' @param aggregationFactor the dataframe of factor to take into account for the randomization
#' @param randomizationNumber the number of random abundance matrix to generate
#' @param seed the seed of the pseudo random number generator
#' @param weightsDataframeRDS the path and name of the file to save generated matrix
#' @param doParallel Indicates if we use parallelism to construct the random matrix
#' @importFrom foreach %dopar%
#' @export

generateRandomMatrix <- function(weights,
                                 aggregationFactor = NULL,
                                 randomizationNumber,
                                 seed = 123456,
                                 weightsDataframeRDS = NULL,
                                 doParallel = TRUE) {
  
  if (!is.null(aggregationFactor) && nrow(weights) != nrow(aggregationFactor)) {
    stop("weights and aggregationFactor must have the same number of rows !")
  }
  
  # Construct the id for aggregation
  if (!is.null(aggregationFactor)) {
    if (is.data.frame(aggregationFactor)) {
      aggregationId <- apply(aggregationFactor, 1, paste, collapse = "_")
    }else {
      aggregationId <- aggregationFactor
    }
  }else {
    # not working with empty id
    aggregationId <- rep(x = "_", times = nrow(weights))
  }
  
  # Construct a list which contains for each aggregation factor the valid weight index,
  # i.e. the sum of weight is not equal to 0
  aggregationFactorIndexList <- list()
  for (agFactor in unique(aggregationId)) {
    aggregationFactorIndexList[[agFactor]] <-
      as.vector(which(colSums(weights[which(aggregationId == agFactor), ]) != 0))
  }
  
  # Set seed for the Pseudo Random Number Generator
  set.seed(seed = seed)
  
  # Initialization of the parallelization if doParallel is true
  if (doParallel == TRUE) {
    nc <- parallel::detectCores()
    cl <- parallel::makeCluster(nc)
    on.exit(expr = parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
  }
  
  # For lintr fake declaration
  randNumber <- 0
  
  # Generation of the random matrix
  weightsDataframe <- foreach::foreach(randNumber = 0:randomizationNumber, .combine = "rbind") %dopar% {
    # Creation of the dataframe which receive the random weights (regarding aggregation factor)
    dataframeToReturn <- data.frame(matrix(data = 0, nrow = nrow(weights), ncol = ncol(weights) + 1))
    colnames(dataframeToReturn) <- c("Number", paste0("Index", seq_len(ncol(weights))))
    dataframeToReturn$Number <- randNumber
    
    # if randNumber = 0, put the original weights data,
    # otherwise weights are shuffle randomly regarding valid index
    if (randNumber == 0) {
      dataframeToReturn[seq_len(nrow(weights)), 2:ncol(dataframeToReturn)] <- weights
    }else {
      for (weightsLineNumber in seq_len(nrow(weights))){
        index <- aggregationFactorIndexList[[aggregationId[weightsLineNumber]]]
        dataframeToReturn[weightsLineNumber, 1 + index] <-
          weights[weightsLineNumber, sample(index, replace = FALSE)]
      }
    }
    return(dataframeToReturn)
  }
  # save the result
  if (!is.null(weightsDataframeRDS)) {
    saveRDS(weightsDataframe, file = weightsDataframeRDS)
  }
  
  return(weightsDataframe)
}

### b. Compare a value to random values ----

#' @title Compare a value to random values
#' @description Compute different statistics (standardized by the distribution of random  values).
#' @concept Statistics
#' @param observedValue the observed value
#' @param randomValues the random Values
#' @param significanceThreshold the array of values used to compute the quantile (c(0.025, 0.975) by default)
#' @return a list corresponding to :
#' - the observed value
#' - quantile values (minimum significance threshold)
#' - quantile values (maximum significance threshold)
#' - significance (observed value not in quantile values)
#' @examples nullModelDistributionStatistics(observedValue = 2, randomValues = c(1, 4, 5, 6, 8),
#' significanceThreshold = c(0.025,0.975))
#' @export

nullModelDistributionStatistics <- function(observedValue,
                                            randomValues,
                                            significanceThreshold = c(0.05, 0.95)) {
  meanRandom <- mean(randomValues, na.rm = T)
  sdRandom <- stats::sd(randomValues, na.rm = T)
  standardizedObserved  <- (observedValue - meanRandom) / sdRandom
  quant <- stats::quantile(x = randomValues,
                           probs = significanceThreshold, 
                           na.rm = T)
  return(list(standardizedObserved,
              (quant[[1]] - meanRandom) / sdRandom,
              (quant[[2]] - meanRandom) / sdRandom,
              observedValue > quant[[2]] || observedValue < quant[[1]]))
}

## Step 2. Function: SKR ANALYSIS ----

#' @title Launch the analysis of the TADs
#' @description Launch the SKR analysis of the TADs, and generate output dataset
#' @concept tad
#' @param weights the dataframe of weights, one row correspond to a series of observation
#' @param weightsFactor the dataframe which contains the different factor linked to the weights
#' @param dataToTreat a vector of the data linked to the different factor
#' @param randomizationNumber the number of random abundance matrix to generate
#' @param aggregationFactorName vector of factor name for the generation of random matrix
#' @param statisticsFactorName vector of factor name for the computation of statistics for each generated matrix
#' @param seed the seed of the pseudo random number generator
#' @param abundanceDataFrameRDS the path and name of the RDS file to load/save the dataframe which
#' contains the observed data and the generated matrix
#' @param weightedMomentsDataFrameRDS the path and name of the RDS file to load/save the dataframe which
#' contains the calculated moments
#' @param statPerObsDataFrameRDS the path and name of the RDS file to load/save the dataframe which
#' contains the statistics for each observed row regarding the random ones
#' @param statPerRandDataFrameRDS the path and name of the RDS file to load/save the dataframe which
#' contains the statistics for each random matrix generated
#' @param statSKRparam the path and name of the RDS file to load/save the dataframe which
#' contains the SKR statistics
#' @param regenerateAbundanceDataFrame boolean to specify if the abundance dataframe is computed again
#' @param regenerateWeightedMomentsDataFrame boolean to specify if the weighted moments dataframe is computed again
#' @param regenerateStatPerObsDataFrame boolean to specify if
#' the statistics per observation dataframe is computed again
#' @param regenerateStatPerRandDataFrame boolean to specify if
#' the statistics per random matrix dataframe is computed again
#' @param significativityThreshold the significance threshold to
#' consider that the observed value is in the randomized value
#' @param doParallel Indicates if we use parallelism to construct the random matrix
#' @param slope_distance Indicates if we use parallelism to construct the random matrix
#' @param intercept_distance Indicates if we use parallelism to construct the random matrix
#' @param lin_mod Indicates the type of linear model to use for (SKR): choose "lm" or "mblm"
#' @export

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
    significativityThreshold = c(0.05, 0.95),
    doParallel = TRUE,
    lin_mod = "lm",
    slope_distance = 1,
    intercept_distance = 1.86
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
    weightedMoments$distanceLaw<- weightedMoments$kurtosis - (slope_distance*weightedMoments$skewness*weightedMoments$skewness + intercept_distance)
    
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
        nullModelDistributionStatistics(
          observedValue = weightedMoments$mean[i],
          randomValues = weightedMoments$mean[(1:randomizationNumber) * nrow(statisticsPerObservation) + i],
          significanceThreshold = significativityThreshold)
      statisticsPerObservation[i, (ncol(weightsFactor) + 5):(ncol(weightsFactor) + 8)] <-
        nullModelDistributionStatistics(
          observedValue = weightedMoments$variance[i],
          randomValues = weightedMoments$variance[(1:randomizationNumber) * nrow(statisticsPerObservation) + i],
          significanceThreshold = significativityThreshold)
      statisticsPerObservation[i, (ncol(weightsFactor) + 9):(ncol(weightsFactor) + 12)] <-
        nullModelDistributionStatistics(
          observedValue = weightedMoments$skewness[i],
          randomValues = weightedMoments$skewness[(1:randomizationNumber) * nrow(statisticsPerObservation) + i],
          significanceThreshold = significativityThreshold)
      statisticsPerObservation[i, (ncol(weightsFactor) + 13):(ncol(weightsFactor) + 16)] <-
        nullModelDistributionStatistics(
          observedValue = weightedMoments$kurtosis[i],
          randomValues = weightedMoments$kurtosis[(1:randomizationNumber) * nrow(statisticsPerObservation) + i],
          significanceThreshold = significativityThreshold)
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
    abundanceDataframe$distanceLaw <- weightedMoments$distanceLaw
    
    for (i in 0:randomizationNumber) {
      for (j in 1:lengthFactor) {
        statisticsPerRandom$Number[i * lengthFactor + j] <- i
        
        statisticsPerRandom[i * lengthFactor + j, statisticsFactorName] <-
          weightsFactor[which(statisticsId == names(statisticsFactorSpeciesList)[j])[1], statisticsFactorName]
        
        dfToAnalyze <- abundanceDataframe[which(x = abundanceDataframe$Number == i), ]
        dfToAnalyze <- dfToAnalyze[which(x = statisticsId == names(statisticsFactorSpeciesList)[j]), ]
        y <- dfToAnalyze$kurtosis
        x <- dfToAnalyze$skewness^2
        distLaw <- dfToAnalyze$distanceLaw^2
        
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
        statisticsPerRandom$RMSE[i * lengthFactor + j] <- sqrt(mean(fit$residuals^2, na.rm = TRUE))
        statisticsPerRandom$Mean_distLaw[i * lengthFactor + j] <- sqrt(mean(distLaw, na.rm = T))
        statisticsPerRandom$CV_distLaw[i * lengthFactor + j] <- sd(distLaw, na.rm = T)*100/mean(distLaw, na.rm = T)
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
                       Slope_SES = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Slope,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Slope,
                         significanceThreshold = significativityThreshold)[[1]][1],
                       Slope_Signi = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Slope,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Slope,
                         significanceThreshold = significativityThreshold)[[4]][1],
                       Intercept_SES = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Intercept,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Intercept,
                         significanceThreshold = significativityThreshold)[[1]][1],
                       Intercept_Signi = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Intercept,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Intercept,
                         significanceThreshold = significativityThreshold)[[4]][1],
                       Rsquare_SES = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Rsquare,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Rsquare,
                         significanceThreshold = significativityThreshold)[[1]][1],
                       Rsquare_Signi = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Rsquare,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Rsquare,
                         significanceThreshold = significativityThreshold)[[4]][1],
                       RMSE_SES = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$RMSE,
                         randomValues = SKRparam[SKRparam$Number > 0,]$RMSE,
                         significanceThreshold = significativityThreshold)[[1]][1],
                       RMSE_Signi = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$RMSE,
                         randomValues = SKRparam[SKRparam$Number > 0,]$RMSE,
                         significanceThreshold = significativityThreshold)[[4]][1],
                       Mean_distLaw_SES = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Mean_distLaw,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Mean_distLaw,
                         significanceThreshold = significativityThreshold)[[1]][1],
                       Mean_distLaw_Signi = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$Mean_distLaw,
                         randomValues = SKRparam[SKRparam$Number > 0,]$Mean_distLaw,
                         significanceThreshold = significativityThreshold)[[4]][1],
                       CV_distLaw_SES = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$CV_distLaw,
                         randomValues = SKRparam[SKRparam$Number > 0,]$CV_distLaw,
                         significanceThreshold = significativityThreshold)[[1]][1],
                       CV_distLaw_Signi = nullModelDistributionStatistics(
                         observedValue = SKRparam[SKRparam[[statisticsFactorName]] == i & SKRparam$Number == 0,]$CV_distLaw,
                         randomValues = SKRparam[SKRparam$Number > 0,]$CV_distLaw,
                         significanceThreshold = significativityThreshold)[[4]][1],
                       statisticsFactorName = i))
  }
  names(SES_SKR)[names(SES_SKR) == "statisticsFactorName"] <- statisticsFactorName
  
  saveRDS(object = SES_SKR, file = statSKRparam)
}

## Step 3. Graphical representation functions ----
### a. Distribution moments representation (mean, variance, skewness and kurtosis) ----

#' @title Graph: Distribution moments
#' @param MOM Moments data frame (mean, variance, skewness, kurtosis)
#' @param SESMOM SES of the Moments data frame and significance compared to null model
#' @param statisticsFactorName column of data use for colors discrimination
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param saveGraphMoments The path to save the graph
#' @export

GraphMoments <- function(MOM,
                         SESMOM,
                         statisticsFactorName,
                         statisticsFactorNameBreaks = NULL,
                         statisticsFactorNameCol = palette(),
                         saveGraphMoments) {
  ggplot2::ggsave(saveGraphMoments,
                  ggpubr::ggarrange(
                    ggplot2::ggplot()+
                      ggplot2::geom_boxplot(data = MOM %>%
                                              dplyr::filter(Number > 0), 
                                            ggplot2::aes(x = "Mean", y = mean),
                                            col = "black", fill = "lightgrey", alpha = 0.4)+
                      ggplot2::geom_point(data = MOM %>% 
                                            dplyr::filter(Number == 0), 
                                          ggplot2::aes(x = "Mean", y = mean, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = "Moments")+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_text(size = 30),
                                     axis.text.x = ggplot2::element_blank(),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"), 
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_boxplot(data = MOM %>% 
                                              dplyr::filter(Number > 0), 
                                            ggplot2::aes(x = "Variance", y = variance),
                                            col = "black", fill = "lightgrey", alpha = 0.4)+
                      ggplot2::geom_point(data = MOM %>% 
                                            dplyr::filter(Number == 0), 
                                          ggplot2::aes(x = "Variance", y = variance, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = "Moments")+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_blank(),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"), 
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_boxplot(data = MOM %>% 
                                              dplyr::filter(Number > 0), 
                                            ggplot2::aes(x = "Skewness", y = skewness),
                                            col = "black", fill = "lightgrey", alpha = 0.4)+
                      ggplot2::geom_point(data = MOM %>% 
                                            dplyr::filter(Number == 0), 
                                          ggplot2::aes(x = "Skewness", y = skewness, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = "Moments")+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_blank(),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"), 
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_boxplot(data = MOM %>% 
                                              dplyr::filter(Number > 0), 
                                            ggplot2::aes(x = "Kurtosis", y = kurtosis),
                                            col = "black", fill = "lightgrey", alpha = 0.4)+
                      ggplot2::geom_point(data = MOM %>% 
                                            dplyr::filter(Number == 0), 
                                          ggplot2::aes(x = "Kurtosis", y = kurtosis, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.4, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = "Moments")+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_blank(),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"), 
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
                      ggplot2::geom_point(data = SESMOM %>% 
                                            dplyr::filter(significanceMean == "TRUE"), 
                                          ggplot2::aes(x = "Mean", y = standardizedObservedMean, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
                      ggplot2::geom_point(data = SESMOM %>% 
                                            dplyr::filter(significanceMean == "FALSE"), 
                                          ggplot2::aes(x = "Mean", y = standardizedObservedMean, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs(y = paste0("SES ", "Moments"))+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_text(size = 30),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"), 
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
                      ggplot2::geom_point(data = SESMOM %>% 
                                            dplyr::filter(significanceVariance == "TRUE"), 
                                          ggplot2::aes(x = "Variance", y = standardizedObservedVariance, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
                      ggplot2::geom_point(data = SESMOM %>% 
                                            dplyr::filter(significanceVariance == "FALSE"), 
                                          ggplot2::aes(x = "Variance", y = standardizedObservedVariance, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs()+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"), 
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
                      ggplot2::geom_point(data = SESMOM %>% 
                                            dplyr::filter(significanceSkewness == "TRUE"), 
                                          ggplot2::aes(x = "Skewness", y = standardizedObservedSkewness, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
                      ggplot2::geom_point(data = SESMOM %>% 
                                            dplyr::filter(significanceSkewness == "FALSE"), 
                                          ggplot2::aes(x = "Skewness", y = standardizedObservedSkewness, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs()+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"), 
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ggplot2::ggplot()+
                      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
                      ggplot2::geom_point(data = SESMOM %>% 
                                            dplyr::filter(significanceKurtosis == "TRUE"), 
                                          ggplot2::aes(x = "Kurtosis", y = standardizedObservedKurtosis, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.8, position = "jitter")+
                      ggplot2::geom_point(data = SESMOM %>% 
                                            dplyr::filter(significanceKurtosis == "FALSE"), 
                                          ggplot2::aes(x = "Kurtosis", y = standardizedObservedKurtosis, col = !!sym(statisticsFactorName), fill = !!sym(statisticsFactorName)), 
                                          shape = 21, size = 4, alpha = 0.2, position = "jitter")+
                      ggplot2::scale_fill_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::scale_color_manual(values = statisticsFactorNameCol, limits = statisticsFactorNameBreaks)+
                      ggplot2::theme_bw()+
                      ggplot2::labs()+
                      ggplot2::theme(plot.title = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_text(size = 20),
                                     axis.title.y = ggplot2::element_blank(),
                                     axis.text.x = ggplot2::element_text(size = 20),
                                     axis.title.x = ggplot2::element_blank(),
                                     legend.title = ggplot2::element_text(size = 30, face = "bold"), 
                                     legend.text = ggplot2::element_text(size = 30),
                                     legend.key.size = ggplot2::unit(1.5, 'cm')),
                    ncol = 4,
                    nrow = 2,
                    common.legend = T,
                    legend = "bottom"
                  ),
                  dpi = 600,
                  width = 15,
                  height = 8
  )
}

### b. SKR representation ----

#' @title Graph: SKR
#' @param MOM moments data frame (mean, variance, skewness, kurtosis) 
#' @param statisticsFactorName column of data use for colors discrimination
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param slope_distance slope of the theoretical distribution law (default: slope = 1 intercept = 1.86 skew-uniform)
#' @param intercept_distance intercept of the theoretical distribution law (default: slope = 1 intercept = 1.86 skew-uniform)
#' @param saveGraphSKR The path to save the graph
#' @export

GraphSKR <- function(
    MOM,
    statisticsFactorName,
    statisticsFactorNameBreaks = NULL,
    statisticsFactorNameCol = palette(),
    slope_distance = 1,
    intercept_distance = 1.86,
    saveGraphSKR
) {
  ggplot2::ggsave(saveGraphSKR,
                  ggplot2::ggplot() +
                    ggplot2::geom_point(data = MOM %>% 
                                          dplyr::filter(Number > 0),
                                        ggplot2::aes(x = skewness**2, y = kurtosis), 
                                        shape = 21, size = 2, alpha = 0.4, col = "#D3D3D3", fill = "#D3D3D3")+
                    ggplot2::geom_smooth(data = MOM %>% 
                                           dplyr::filter(Number > 0),
                                         ggplot2::aes(x = skewness**2, y = kurtosis, group = Number), 
                                         col = "#D3D3D3", fill = "#D3D3D3", se = F, method = "lm", formula = y ~ x, linetype = 1, linewidth = 0.5, alpha = 0.1)+
                    ggplot2::geom_abline(intercept = intercept_distance, slope = slope_distance, linetype = "dashed", linewidth = 2) +
                    ggplot2::geom_point(data = MOM %>% 
                                          dplyr::filter(Number == 0),
                                        ggplot2::aes(x = skewness**2, y = kurtosis, fill = !!rlang::sym(statisticsFactorName)), 
                                        shape = 21, size = 6, alpha = 0.4)+
                    ggplot2::geom_smooth(data = MOM %>% 
                                           dplyr::filter(Number == 0),
                                         ggplot2::aes(x = skewness**2, y = kurtosis, col = !!rlang::sym(statisticsFactorName), fill = !!rlang::sym(statisticsFactorName)), 
                                         se = F, method = "lm", formula = y ~ x, linetype = 1, linewidth = 2, alpha = 0.1)+
                    ggpubr::stat_regline_equation(data = MOM %>% 
                                                    dplyr::filter(Number == 0),
                                                  ggplot2::aes(x = skewness**2, y = kurtosis, col = !!rlang::sym(statisticsFactorName)),
                                                  alpha = 1, size = 8)+
                    ggplot2::scale_fill_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
                    ggplot2::scale_color_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
                    ggplot2::xlim(0, 10)+
                    ggplot2::ylim(0, 20)+
                    ggplot2::theme_bw()+
                    ggplot2::labs(x = "Skewness²", y = "Kurtosis")+
                    ggplot2::theme(legend.position = "bottom",
                                   plot.title = ggplot2::element_blank(),
                                   axis.text.y = ggplot2::element_text(size = 35),
                                   axis.title.y = ggplot2::element_text(size = 40),
                                   axis.title.x = ggplot2::element_text(size = 40),
                                   axis.text.x = ggplot2::element_text(size = 35)),
                  dpi = 600,
                  height = 10,
                  width = 10
  )
}

### c. Parameters of the SKR ----

#' @title Graph: parameters of the SKR
#' @param SKRparam SES of SKR parameters data frame (SES and Significance)
#' @param statisticsFactorName column of data use for colors discrimination
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param slope_distance slope of the theoretical distribution law (default: slope = 1 intercept = 1.86 skew-uniform)
#' @param intercept_distance intercept of the theoretical distribution law (default: slope = 1 intercept = 1.86 skew-uniform)
#' @param saveGraphparamSKR The path to save the graph
#' @export

GraphparamSKR <- function(SKRparam,
                          statisticsFactorName,
                          statisticsFactorNameBreaks = NULL,
                          statisticsFactorNameCol = palette(),
                          slope_distance = 1,
                          intercept_distance = 1.86,
                          saveGraphparamSKR) {
  
  if(slope_distance == 1 & intercept_distance == 1.86){
    title_dist_law <- "distance to 
Skew-Uniform"
  }else if (slope_distance == 1 & intercept_distance == 1){
    title_dist_law <- "distance to 
Lower Boundary"
  }else{
    title_dist_law <- paste0("distance to 
K = ", slope_distance, " x S² + ", intercept_distance)
  }
  
  ggplot2::ggsave(
    saveGraphparamSKR,
    ggplot2::ggplot()+
      ggplot2::geom_abline(intercept = 0, slope = 0, color = "grey", linewidth = 1, linetype = "dashed")+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(Slope_Signi == TRUE), 
                          ggplot2::aes(x =  "Slope", y = Slope_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(Slope_Signi == FALSE), 
                          ggplot2::aes(x =  "Slope", y = Slope_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(Intercept_Signi == TRUE), 
                          ggplot2::aes(x =  "Intercept", y = Intercept_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(Intercept_Signi == FALSE), 
                          ggplot2::aes(x =  "Intercept", y = Intercept_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(Rsquare_Signi == TRUE), 
                          ggplot2::aes(x =  "R²", y = Rsquare_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(Rsquare_Signi == FALSE), 
                          ggplot2::aes(x =  "R²", y = Rsquare_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(Mean_distLaw_Signi == TRUE), 
                          ggplot2::aes(x = title_dist_law, y = Mean_distLaw_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(Mean_distLaw_Signi == FALSE), 
                          ggplot2::aes(x = title_dist_law, y = Mean_distLaw_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(CV_distLaw_Signi == TRUE), 
                          ggplot2::aes(x = paste0("CV ", title_dist_law), y = CV_distLaw_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(CV_distLaw_Signi == FALSE), 
                          ggplot2::aes(x = paste0("CV ", title_dist_law), y = CV_distLaw_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(RMSE_Signi == TRUE), 
                          ggplot2::aes(x =  "RMSE", y = RMSE_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.8, size = 6, color = "black", shape = 21)+
      ggplot2::geom_point(data = SKRparam %>% 
                            dplyr::filter(RMSE_Signi == FALSE), 
                          ggplot2::aes(x =  "RMSE", y = RMSE_SES, fill = !!sym(statisticsFactorName)), 
                          alpha = 0.2, size = 6, color = "black", shape = 21)+
      ggplot2::scale_x_discrete(limits = c("Slope", "Intercept", "R²", "RMSE", title_dist_law, paste0("CV ", title_dist_law)))+
      ggplot2::scale_fill_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
      ggplot2::scale_color_manual(limits = statisticsFactorNameBreaks, values = statisticsFactorNameCol)+
      ggplot2::theme_bw()+
      ggplot2::labs(title = paste0("Parameters of the SKR"), y = "SES")+
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
                     axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = 10, face = "bold")),
    dpi = 600,
    width = 15,
    height = 5
  )
}
