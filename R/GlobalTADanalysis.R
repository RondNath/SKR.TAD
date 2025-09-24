## Step 4. Global TADs analysis ----

#' @title Global laucher for TADs analysis (group all functions of the package)
#' @description Launch the SKR analysis of the TADs, generate output dataset (distribution moments, SKR parameters, and statistics) and related graphics
#' @param Abundance the dataframe of abundance (or related weights measure), one row correspond to a series of observation
#' @param AbundanceFactor the dataframe which contains the different factor linked to the abundance dataframe
#' @param TraitData a vector of the trait data linked to the different abundance factor
#' @param randomizationNumber the number of random abundance matrix to generate
#' @param randomizationFactorName vector of factor name for the generation of random matrix
#' @param statisticsFactorName vector of factor name for the computation of statistics for each generated matrix
#' @param statisticsFactorNameBreaks vector of factor levels of the statisticsFactorName, same dimension than statisticsFactorNameCol
#' @param statisticsFactorNameCol vector of colors, same dimension than statisticsFactorNameBreaks
#' @param seed the seed of the pseudo random number generator
#' @param path_abundanceDataFrame the path and name of the RDS file to load/save the dataframe which contains the observed abundance data and the generated random abundance matrix
#' @param path_MomentsDataFrame the path and name of the RDS file to load/save the dataframe which contains the calculated moments (for observations and randomizations)
#' @param path_SESMomentsDataFrame the path and name of the RDS file to load/save the dataframe which contains the moments standardized effect size (observations regarding randomizations)
#' @param path_SKRDataFrame the path and name of the RDS file to load/save the dataframe which contains the SKR parameters (for observations and randomizations)
#' @param path_SESSKRDataFrame the path and name of the RDS file to load/save the dataframe which contains the SKR parameters standardized effect size (observations regarding randomizations)
#' @param path_GraphMoments The path to save the graph of the moments (mean, variance, skewness and kurtosis)
#' @param path_GraphSKR The path to save the graph of the SKR
#' @param path_GraphparamSKR The path to save the graph of the SKR parameters
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
#' @returns
#' RDS files with:
#' abundance for observed and randomized communities,
#' moments (mean, variance, skewness & kurtosis) for observed and randomized communities,
#' Standardized Effect Size (SES) values of the moments for observed compared to randomized communities and significance,
#' SKR parameters for observed and randomized communities
#' SES of the SKR parameters (slope, intercept, Rsquare, distance from predicted TADs, distance from reference TADs & CV of the distance from reference TADs).
#' PNG graphics:
#' Plot of the moments in two panels up) the moments, bottom) the SES moments values
#' Plot of the SKR for observed and randomized communities in two panels left) kurtosis ~ skewness, right) kurtosis ~ skewness²
#' Plot of the SKR parameters values and significance
#' @export
#' @examples
#' head(data)
#' head(trait)
#' GlobalTADanalysis(
#' Abundance = abundance[,5:102],
#' AbundanceFactor = abundance[,c("Year", "Plot", "Treatment", "Bloc")],
#' TraitData = log(trait[["SLA"]] + 1),
#' randomizationFactorName = c("Year", "Bloc"),
#' statisticsFactorName = c("Treatment"),
#' statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
#' statisticsFactorNameCol = c("#1A85FF", "#D41159"),
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
#' path_GraphMoments = NULL,
#' path_GraphSKR = NULL,
#' path_GraphparamSKR = NULL,
#' significanceThreshold = c(0.05, 0.95),
#' slope_ref_TADs = 1,
#' intercept_ref_TADs = 1.86,
#' distance_metric = "RMSE",
#' lin_mod = "lm",
#' doParallel = FALSE
#' )

GlobalTADanalysis <- function(
    Abundance,
    AbundanceFactor,
    TraitData,
    randomizationNumber,
    randomizationFactorName = NULL,
    statisticsFactorName = NULL,
    statisticsFactorNameBreaks = NULL,
    statisticsFactorNameCol = grDevices::palette(),
    seed = 123456,
    path_abundanceDataFrame = NULL,
    path_MomentsDataFrame = NULL,
    path_SESMomentsDataFrame = NULL,
    path_SKRDataFrame = NULL,
    path_SESSKRDataFrame = NULL,
    path_GraphMoments = NULL,
    path_GraphSKR = NULL,
    path_GraphparamSKR = NULL,
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
    distance_metric = "RMSE"
){
  return(list(
    DataAnalysisTAD(
      Abundance = Abundance,
      AbundanceFactor = AbundanceFactor,
      TraitData = TraitData,
      randomizationNumber = randomizationNumber,
      randomizationFactorName = randomizationFactorName,
      statisticsFactorName = statisticsFactorName,
      seed = seed,
      path_abundanceDataFrame = path_abundanceDataFrame,
      path_MomentsDataFrame = path_MomentsDataFrame,
      path_SESMomentsDataFrame = path_SESMomentsDataFrame,
      path_SKRDataFrame = path_SKRDataFrame,
      path_SESSKRDataFrame = path_SESSKRDataFrame,
      regenerate_abundanceDataFrame = regenerate_abundanceDataFrame,
      regenerate_MomentsDataFrame = regenerate_MomentsDataFrame,
      regenerate_SESMomentsDataFrame = regenerate_SESMomentsDataFrame,
      regenerate_SKRDataFrame = regenerate_SKRDataFrame,
      regenerate_SESSKRDataFrame = regenerate_SESSKRDataFrame,
      significanceThreshold = significanceThreshold,
      doParallel = doParallel,
      lin_mod = lin_mod,
      slope_ref_TADs = slope_ref_TADs,
      intercept_ref_TADs = intercept_ref_TADs,
      distance_metric = distance_metric
    ),
    GraphMoments(
      MomentsDataFrame = MomentsDataFrame,
      SESMomentsDataFrame = SESMomentsDataFrame,
      statisticsFactorName = statisticsFactorName,
      statisticsFactorNameBreaks = statisticsFactorNameBreaks,
      statisticsFactorNameCol = statisticsFactorNameCol,
      path_GraphMoments = path_GraphMoments
    ),
    GraphSKR(
      MomentsDataFrame = MomentsDataFrame,
      statisticsFactorName = statisticsFactorName,
      statisticsFactorNameBreaks = statisticsFactorNameBreaks,
      statisticsFactorNameCol = statisticsFactorNameCol,
      slope_ref_TADs = slope_ref_TADs,
      intercept_ref_TADs = intercept_ref_TADs,
      path_GraphSKR = path_GraphSKR
    ),
    GraphparamSKR(
      SESSKRDataFrame = SESSKRDataFrame,
      statisticsFactorName = statisticsFactorName,
      statisticsFactorNameBreaks = statisticsFactorNameBreaks,
      statisticsFactorNameCol = statisticsFactorNameCol,
      slope_ref_TADs = slope_ref_TADs,
      intercept_ref_TADs = intercept_ref_TADs,
      path_GraphparamSKR = path_GraphparamSKR
    )))
}
