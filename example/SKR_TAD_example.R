# I. PACKAGES ----
library(devtools)
devtools::install_github("RondNath/SKR.TAD", force = T)
library(SKR.TAD)

# Set working directory to source .R file location used to launch analysis with SKR.TAD package
# "Session" > "Set Working Directory" > "To Source File Location"

# II. Randomizations ----

SKR.TAD::AbundanceRandomization(
  Abundance = SKR.TAD::abundance[,5:102],
  randomizationFactor =  SKR.TAD::abundance[,c("Year", "Bloc")],
  randomizationNumber = 1000,
  seed = 666,
  path_abundanceDataFrame = "./Output/abundanceDataFrame.RDS",
  doParallel = TRUE
)

# III. LAUNCH SKR ANALYSIS ----

SKR.TAD::DataAnalysisTAD(Abundance = SKR.TAD::abundance[,5:102],
                         AbundanceFactor = SKR.TAD::abundance[,c("Year", "Plot", "Treatment", "Bloc")],
                         TraitData = log(SKR.TAD::trait[["SLA"]] + 1),
                         randomizationFactorName = c("Year", "Bloc"),
                         statisticsFactorName = c("Treatment"),
                         regenerate_abundanceDataFrame = T,
                         regenerate_MomentsDataFrame = T,
                         regenerate_SESMomentsDataFrame = T,
                         regenerate_SKRDataFrame = T,
                         regenerate_SESSKRDataFrame = T,
                         randomizationNumber = 1000,
                         seed = 666,
                         path_abundanceDataFrame = "./Output/abundanceDataFrame.RDS",
                         path_MomentsDataFrame = "./Output/MomentsDataFrame.RDS",
                         path_SESMomentsDataFrame = "./Output/SESMomentsDataFrame.RDS",
                         path_SKRDataFrame = "./Output/SKRDataFrame.RDS",
                         path_SESSKRDataFrame = "./Output/SESSKRDataFrame.RDS",
                         significanceThreshold = c(0.05, 0.95),
                         slope_ref_TADs = 1,
                         intercept_ref_TADs = 1.86,
                         distance_metric = "RMSE",
                         lin_mod = "lm",
                         doParallel = T)

# III. PLOT SKR RESULTS ----
SKR.TAD::GraphMoments(
  MomentsDataFrame = SKR.TAD::MomentsDataFrame,
  SESMomentsDataFrame = SKR.TAD::SESMomentsDataFrame,
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  path_GraphMoments = "./Output/Moments.png"
)

SKR.TAD::GraphSKR(
  MomentsDataFrame = SKR.TAD::MomentsDataFrame,
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  slope_ref_TADs = 1,
  intercept_ref_TADs = 1.86,
  path_GraphSKR = "./Output/SKR.png"
)

SKR.TAD::GraphparamSKR(
  SESSKRDataFrame = SKR.TAD::SESSKRDataFrame,
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  slope_ref_TADs = 1,
  intercept_ref_TADs = 1.86,
  path_GraphparamSKR = "./Output/SKRparam.png"
)

# IV. GLOBAL TADs ANALYSIS (Group all functions) ----

SKR.TAD::GlobalTADLauncher(
  Abundance = SKR.TAD::abundance[,5:102],
  AbundanceFactor = SKR.TAD::abundance[,c("Year", "Plot", "Treatment", "Bloc")],
  TraitData = log(SKR.TAD::trait[["SLA"]] + 1),
  randomizationFactorName = c("Year", "Bloc"),
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  regenerate_abundanceDataFrame = T,
  regenerate_MomentsDataFrame = T,
  regenerate_SESMomentsDataFrame = T,
  regenerate_SKRDataFrame = T,
  regenerate_SESSKRDataFrame = T,
  randomizationNumber = 1000,
  seed = 666,
  path_abundanceDataFrame = "./Output/abundanceDataFrame.RDS",
  path_MomentsDataFrame = "./Output/MomentsDataFrame.RDS",
  path_SESMomentsDataFrame = "./Output/SESMomentsDataFrame.RDS",
  path_SKRDataFrame = "./Output/SKRDataFrame.RDS",
  path_SESSKRDataFrame = "./Output/SESSKRDataFrame.RDS",
  path_GraphMoments = "./Output/Moments.png",
  path_GraphSKR = "./Output/SKR.png",
  path_GraphparamSKR = "./Output/SKRparam.png",
  significanceThreshold = c(0.05, 0.95),
  slope_ref_TADs = 1,
  intercept_ref_TADs = 1.86,
  distance_metric = "RMSE",
  lin_mod = "lm",
  doParallel = T
)
