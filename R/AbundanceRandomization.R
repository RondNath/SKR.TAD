## Step 1. FUNCTION: RANDOMIZATION ----
### a. Abundance Randomization ----

#' @title Abundance Randomization
#' @description Generate and save random matrix (abundance randomization)
#' @param Abundance the dataframe of abundance (or related weights measure), one row correspond to a series of observation
#' @param randomizationFactor vector of factor name for the generation of random matrix
#' @param randomizationNumber the number of random abundance matrix to generate
#' @param seed the seed of the pseudo random number generator
#' @param path_abundanceDataFrame the path and name of the RDS file to load/save the dataframe which contains the observed abundance data and the generated random abundance matrix
#' @param doParallel Indicates if we use parallelism to construct the random matrix
#' @importFrom foreach %dopar%
#' @returns RDS file with the abundance of observed and randomized communities
#' @export
#' @examples
#' head(abundance)
#' AbundanceRandomization(
#' Abundance = abundance[,5:102],
#' randomizationFactor = abundance[,c("Year", "Bloc")],
#' randomizationNumber = 1000,
#' seed = 666,
#' path_abundanceDataFrame = NULL,
#' doParallel = FALSE
#' )

AbundanceRandomization <- function(
    Abundance,
    randomizationFactor = NULL,
    randomizationNumber,
    seed = 123456,
    path_abundanceDataFrame = NULL,
    doParallel = TRUE
) {

  if (!is.null(randomizationFactor) && nrow(Abundance) != nrow(randomizationFactor)) {
    stop("Abundance and randomizationFactor must have the same number of rows !")
  }

  # Construct the id for aggregation
  if (!is.null(randomizationFactor)) {
    if (is.data.frame(randomizationFactor)) {
      aggregationId <- apply(randomizationFactor, 1, paste, collapse = "_")
    }else {
      aggregationId <- randomizationFactor
    }
  }else {
    # not working with empty id
    aggregationId <- rep(x = "_", times = nrow(Abundance))
  }

  # Construct a list which contains for each randomization factor the valid weight (abundance) index,
  # i.e. the sum of weight is not equal to 0
  randomizationFactorIndexList <- list()
  for (agFactor in unique(aggregationId)) {
    randomizationFactorIndexList[[agFactor]] <-
      as.vector(which(colSums(Abundance[which(aggregationId == agFactor), ]) != 0))
  }

  # Set seed for the Pseudo Random Number Generator
  set.seed(seed = seed)

  # Initialization of the parallelization if doParallel is true
  if (doParallel == TRUE) {
    nc <- parallel::detectCores()
    cl <- parallel::makeCluster(nc)
    on.exit(expr = parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
  }else{
    nc <- 1
    cl <- parallel::makeCluster(nc)
    on.exit(expr = parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
  }

  # For lintr fake declaration
  randNumber <- 0

  # Generation of the random matrix
  abundanceDataframe <- foreach::foreach(randNumber = 0:randomizationNumber, .combine = "rbind") %dopar% {
    # Creation of the dataframe which receive the random weights (regarding aggregation factor)
    dataframeToReturn <- data.frame(matrix(data = 0, nrow = nrow(Abundance), ncol = ncol(Abundance) + 1))
    colnames(dataframeToReturn) <- c("Number", paste0("Index", seq_len(ncol(Abundance))))
    dataframeToReturn$Number <- randNumber

    # if randNumber = 0, put the original weights data,
    # otherwise weights are shuffle randomly regarding valid index
    if (randNumber == 0) {
      dataframeToReturn[seq_len(nrow(Abundance)), 2:ncol(dataframeToReturn)] <- Abundance
    }else {
      for (weightsLineNumber in seq_len(nrow(Abundance))){
        index <- randomizationFactorIndexList[[aggregationId[weightsLineNumber]]]
        dataframeToReturn[weightsLineNumber, 1 + index] <-
          Abundance[weightsLineNumber, sample(index, replace = FALSE)]
      }
    }
    return(dataframeToReturn)
  }
  # save the result
  if (!is.null(path_abundanceDataFrame)) {
    saveRDS(abundanceDataframe, file = path_abundanceDataFrame)
  }

  return(abundanceDataframe)
}
