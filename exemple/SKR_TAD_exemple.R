# I. PACKAGES ----
library(devtools)
devtools::install_github("RondNath/SKR.TAD", force = T)
library(SKR.TAD)

# II. LAUNCH SKR ANALYSIS ----

# Load abundance Data
abundance <- r4urep::loadTableFile(path = "./Input/abundance.csv",
                            colFactor = c("Plot", "Treatment", "Year", "Bloc"),
                            sep = ";",
                            dec = ".")
# Load trait Data
trait <- r4urep::loadTableFile(path = "./Input/trait.csv",
                               sep = ";",
                               dec= ".")

SKR.TAD::DataAnalysisTAD(
  weights = abundance[,5:102],
  weightsFactor = abundance[,c("Year", "Plot", "Treatment", "Bloc")],
  dataToTreat = log(trait[["SLA"]]),
  aggregationFactorName = c("Year", "Bloc"),
  statisticsFactorName = c("Treatment"),
  regenerateAbundanceDataFrame = T,
  regenerateWeightedMomentsDataFrame = T,
  regenerateStatPerObsDataFrame = T,
  regenerateStatPerRandDataFrame = T,
  randomizationNumber = 1000,
  seed = 666,
  abundanceDataFrameRDS = "./Output/abundanceDataFrame.RDS",
  weightedMomentsDataFrameRDS = "./Output/MomentsDataFrame.RDS",
  statPerObsDataFrameRDS = "./Output/SES_MomentsDataFrame.RDS",
  statPerRandDataFrameRDS = "./Output/SKRDataFrame.RDS",
  statSKRparam = "./Output/SES_SKRDataFrame.RDS",
  significanceThreshold = c(0.05, 0.95),
  slope_speTADs = 1,
  intercept_speTADs = 1.86,
  distance_metric = "RMSE",
  lin_mod = "lm"
)

# III. PLOT SKR RESULTS ----

SKR.TAD::GraphMoments(
  MOM = readRDS("./Output/MomentsDataFrame.RDS"),
  SESMOM = readRDS("./Output/SES_MomentsDataFrame.RDS"),
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  saveGraphMoments = "./Output/Moments.png"
)

SKR.TAD::GraphSKR(
  MOM = readRDS("./Output/MomentsDataFrame.RDS"),
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  slope_speTADs = 1,
  intercept_speTADs = 1.86,
  saveGraphSKR = "./Output/SKR.png"
)

SKR.TAD::GraphparamSKR(
  SKRparam = readRDS("./Output/SES_SKRDataFrame.RDS"),
  statisticsFactorName = c("Treatment"),
  statisticsFactorNameBreaks = c("Mown_Unfertilized", "Mown_NPK"),
  statisticsFactorNameCol = c("#1A85FF", "#D41159"),
  slope_speTADs = 1,
  intercept_speTADs = 1.86,
  saveGraphparamSKR = "./Output/paramSKR.png"
)
