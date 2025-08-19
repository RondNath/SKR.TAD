## Step 1. FUNCTION: RANDOMIZATION ----
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
