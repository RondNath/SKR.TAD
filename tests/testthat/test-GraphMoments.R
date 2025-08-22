library(testthat)
library(ggplot2)
library(dplyr)
library(rlang)
library(ggpubr)
library(SKR.TAD)

MOM = tibble(
  Number = c(0, 0, 1, 1),
  Factor = c("A", "A", "B", "B"),
  mean = c(1, 2, 3, 4),
  variance = c(0.5, 1.5, 2.5, 3.5),
  skewness = c(0, 0.1, 0.2, 0.3),
  kurtosis = c(3, 3.1, 3.2, 3.3)
)

SESMOM = tibble(
  Factor = c("A", "A", "B", "B"),
  significanceMean = c(TRUE, FALSE, TRUE, FALSE),
  significanceVariance = c(TRUE, FALSE, TRUE, FALSE),
  significanceSkewness = c(TRUE, FALSE, TRUE, FALSE),
  significanceKurtosis = c(TRUE, FALSE, TRUE, FALSE),
  standardizedObservedMean = c(1, 1.5, 2, 2.5),
  standardizedObservedVariance = c(0.2, 0.3, 0.4, 0.5),
  standardizedObservedSkewness = c(0.1, 0.2, 0.3, 0.4),
  standardizedObservedKurtosis = c(0.1, 0.2, 0.3, 0.4)
)

# 1/ Test of graphic saving ----
testthat::test_that("GraphMoments build png graphic", {

  # Create temporary file
  temp_file <- tempfile(fileext = ".png")

  GraphMoments(MOM = MOM,
  SESMOM = SESMOM,
  statisticsFactorName = "Factor",
  statisticsFactorNameBreaks = c("A", "B"),
  statisticsFactorNameCol = c("pink", "cyan"),
  saveGraphMoments = temp_file)

  # Check whether the file has been created
  testthat::expect_true(file.exists(temp_file))

  # Check whether the file is png
  testthat::expect_true(tools::file_ext(temp_file) == "png")

  # Delete temporary file
  unlink(temp_file)
})

# 2/ Test invalid Inputs ----
testthat::test_that("GraphMoments invalid inputs", {
  testthat::expect_error(GraphMoments(MOM = NULL,
                            SESMOM = SESMOM,
                            statisticsFactorName = "Factor",
                            saveGraphMoments = tempfile(fileext = ".png")),
               "argument \"MOM\" is missing")

  testthat::expect_error(GraphMoments(MOM = MOM,
                            SESMOM = NULL,
                            statisticsFactorName = "mean",
                            saveGraphMoments = tempfile(fileext = ".png")),
               "argument \"SESMOM\" is missing")

  testthat::expect_error(GraphMoments(MOM = MOM,
                            SESMOM = SESMOM,
                            statisticsFactorName = NULL,
                            saveGraphMoments = tempfile(fileext = ".png")),
               "argument \"statisticsFactorName\" is missing")
})
