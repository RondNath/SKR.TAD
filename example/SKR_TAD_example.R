# I. PACKAGES ----
library(devtools)
devtools::install_github("RondNath/SKR.TAD", force = T)
library(SKR.TAD)

# Set working directory to source .R file location used to launch analysis with SKR.TAD package
# "Session" > "Set Working Directory" > "To Source File Location"

# II. Randomizations ----

SKR.TAD::AbundanceRandomization(Abundance = SKR.TAD::abundance[,5:102],
  randomizationFactor =  SKR.TAD::abundance[,c("Year", "Bloc")],
  randomizationNumber = 1000,
  seed = 666,
  path_abundanceDataFrame = "./Output/abundanceDataFrame.RDS",
  doParallel = TRUE)

# III. LAUNCH SKR ANALYSIS ----

SKR.TAD::DataAnalysisTAD(Abundance = SKR.TAD::abundance[,5:102],
  AbundanceFactor = SKR.TAD::abundance[,c("Year", "Plot", "Treatment", "Bloc")],
  TraitData = log(SKR.TAD::trait[["SLA"]] + 1),
  randomizationFactorName = c("Year", "Bloc"),
  statisticsFactorName = c("Treatment"),
  regenerate_abundanceDataFrame = T,
  regenerate_momentsDataFrame = T,
  regenerate_SESmomentsDataFrame = T,
  regenerate_SKRDataFrame = T,
  regenerate_SESSKRDataFrame = T,
  randomizationNumber = 100,
  seed = 666,
  path_abundanceDataFrame = "./Output/abundanceDataFrame.RDS",
  path_momentsDataFrame = "./Output/MomentsDataFrame.RDS",
  path_SESmomentsDataFrame = "./Output/SES_MomentsDataFrame.RDS",
  path_SKRDataFrame = "./Output/SKRDataFrame.RDS",
  path_SESSKRDataFrame = "./Output/SES_SKRDataFrame.RDS",
  significanceThreshold = c(0.05, 0.95),
  slope_ref_TADs = 1,
  intercept_ref_TADs = 1.86,
  distance_metric = "RMSE",
  lin_mod = "lm")

# III. PLOT SKR RESULTS ----

SKR.TAD::GraphMoments(moments = readRDS("./Output/MomentsDataFrame.RDS"),
  SESmoments = readRDS("./Output/SES_MomentsDataFrame.RDS"),
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  path_GraphMoments = "./Output/Moments.png")

SKR.TAD::GraphSKR(moments = readRDS("./Output/MomentsDataFrame.RDS"),
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  slope_ref_TADs = 1,
  intercept_ref_TADs = 1.86,
  path_GraphSKR = "./Output/SKR.png")

SKR.TAD::GraphparamSKR(SESSKRparam = readRDS("./Output/SES_SKRDataFrame.RDS"),
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  slope_ref_TADs = 1,
  intercept_ref_TADs = 1.86,
  path_GraphparamSKR = "./Output/paramSKR.png")

# IV. GLOBAL TADs ANALYSIS (Group all functions) ----

SKR.TAD::GlobalTADanalysis(Abundance = SKR.TAD::abundance[,5:102],
                           AbundanceFactor = SKR.TAD::abundance[,c("Year", "Plot", "Treatment", "Bloc")],
                           TraitData = log(SKR.TAD::trait[["SLA"]] + 1),
                           randomizationFactorName = c("Year", "Bloc"),
                           statisticsFactorName = c("Treatment"),
                           statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
                           statisticsFactorNameCol = c("#1A85FF", "#D41159"),
                           regenerate_abundanceDataFrame = T,
                           regenerate_momentsDataFrame = T,
                           regenerate_SESmomentsDataFrame = T,
                           regenerate_SKRDataFrame = T,
                           regenerate_SESSKRDataFrame = T,
                           randomizationNumber = 100,
                           seed = 666,
                           path_abundanceDataFrame = "./Output/abundanceDataFrame.RDS",
                           path_momentsDataFrame = "./Output/MomentsDataFrame.RDS",
                           path_SESmomentsDataFrame = "./Output/SES_MomentsDataFrame.RDS",
                           path_SKRDataFrame = "./Output/SKRDataFrame.RDS",
                           path_SESSKRDataFrame = "./Output/SES_SKRDataFrame.RDS",
                           path_GraphMoments = "./Output/Moments.png",
                           path_GraphSKR = "./Output/SKR.png",
                           path_GraphparamSKR = "./Output/paramSKR.png",
                           significanceThreshold = c(0.05, 0.95),
                           slope_ref_TADs = 1,
                           intercept_ref_TADs = 1.86,
                           distance_metric = "RMSE",
                           lin_mod = "lm")
