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
