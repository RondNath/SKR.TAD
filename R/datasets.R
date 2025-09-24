#' abundance: abundance of observed grasslands plant communities
#'
#' DataFrame of abundance of grassland plant (SP1 to SP110) communities observed over space (Plot) and time (Year) under contrasting management practices (Mown_NPK vs. Mown_Unfertilized) in Massif-central
#'
#' @format DataFrame with 102 columns: `Plot`, `Year`, `Treatment`, `Bloc`, `SP1` to `SP110`
#' @source Observed permanent grassland in Massif-central, ORE long-term experiment (Louault et al. 2017)
"abundance"

#' trait: SLA trait values for grassland species
#'
#' DataFrame of SLA trait values of grassland plant species(SP1 to SP110) observed in the abundance dataset (permanent grassland plant species in Massif-central), trait value from TRY DataBase (Kattge et al. 2020)
#'
#' @format DataFrame with 2 columns: `Species`, `SLA`
#' @source TRY DataBase (Kattge et al. 2020)
"trait"

#' abundanceDataFrame: observed and randomized abundance
#'
#' DataFrame which contains the observed abundance data and the generated random abundance matrix from the initial "abundance" dataframe (observed Massif-central grasslands communities)
#'
#' @format DataFrame with 99 columns: `Number`, `Index1` to `Index98`
#' @source from the function of the package: "AbundanceRandomization"
"abundanceDataFrame"

#' MomentsDataFrame: computed moments (for observations and randomizations)
#'
#' DataFrame which contains the calculated moments (for observations and randomizations), computed from "abundanceDataFrame" (observed and randomized Massif-central grassland communities)
#'
#' @format DataFrame with 10 columns: `Year`, `Plot`, `Treatment`, `Bloc`, `Number`, `mean`, `variance`, `skewness`, `kurtosis`, `distance_reference_TADs`
#' @source from the function of the package: "DataAnalysisTAD"
"MomentsDataFrame"

#' SES_MomentsDataFrame: moments statistics analysis
#'
#' DataFrame which contains the moments standardized effect size (observations regarding randomizations) and statistics, computed from "MomentsDataFrame"
#'
#' @format DataFrame with 20 columns: `Year`, `Plot`, `Treatment`, `Bloc`, `SESMean`, `SESMinQuantileMean`, `SESMaxQuantileMean`, `significanceMean`, `SESVariance`, `SESMinQuantileVariance`, `SESMaxQuantileVariance`, `significanceVariance`, `SESSkewness`, `SESMinQuantileSkewness`, `SESMaxQuantileSkewness`, `significanceSkewness`, `SESKurtosis`, `SESMinQuantileKurtosis`, `SESMaxQuantileKurtosis`, `significanceKurtosis`
#' @source from the function of the package: "DataAnalysisTAD"
"SES_MomentsDataFrame"

#' SKRDataFrame: SKR parameters
#'
#' DataFrame which contains the SKR parameters (for observations and randomizations), computed from "MomentsDataFrame"
#'
#' @format DataFrame with 8 columns: `Number`, `Treatment`, `Slope`, `Intercept`, `Rsquare`, `distance_predicted_TADs`, `distance_reference_TADs`, `CVdistance_reference_TADs`
#' @source from the function of the package: "DataAnalysisTAD"
"SKRDataFrame"

#' SES_SKRDataFrame: SKR parameters statistics analysis
#'
#' DataFrame which contains the SKR parameters standardized effect size (observations regarding randomizations) and statistics, computed from "SKRDataFrame"
#'
#' @format DataFrame with 25 columns: `Treatment`, `SESSlope`, `SESMinQuantileSlope`, `SESMaxQuantileSlope`, `significanceSlope`, `SESIntercept`, `SESMinQuantileIntercept`, `SESMaxQuantileIntercept`, `significanceIntercept`, `SESRsquare`, `SESMinQuantileRsquare`, `SESMaxQuantileRsquare`, `significanceRsquare`, `SESdistance_predicted_TADs`, `SESMinQuantiledistance_predicted_TADs`, `SESMaxQuantiledistance_predicted_TADs`, `significancedistance_predicted_TADs`, `SESdistance_reference_TADs`, `SESMinQuantiledistance_reference_TADs`, `SESMaxQuantiledistance_reference_TADs`, `significancedistance_reference_TADs`, `SESCVdistance_reference_TADs`, `SESMinQuantileCV_distance_reference_TADs`, `SESMaxQuantileCV_distance_reference_TADs`, `significanceCV_distance_reference_TADs`
#' @source from the function of the package: "DataAnalysisTAD"
"SES_SKRDataFrame"
