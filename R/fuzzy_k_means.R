#' Fuzzy k-means clustering
#'
#' @description Performs a vectorized fuzzy k-means clustering, this procedure
#' it is not spatial. The function is almost a wrapper of the function cmeans
#' from the package e1071. Is intended to be used when `KM-sPC` procedure is not
#' possible because data set has only 1 variable.
#'
#' @inheritParams kmspc
#' @return a list with classification results and indices to select best number of
#' clusters.
#' @export
#'
fuzzy_k_means <- function(data,
                          variables,
                          number_cluster = 3:5,
                          fuzzyness = 1.2,
                          distance = "euclidean") {

  if (missing(variables)) {
    myNumVars <-
      unlist(lapply(sf::st_drop_geometry(data), is.numeric))
    if (sum(myNumVars) == 0) {
      stop('Non numeric variables were found in data')
    }
    warning("The numeric Variable will be used to make clusters",
            call. = FALSE)
    variables <- names(sf::st_drop_geometry(data))[myNumVars]
  }

  # if (!inherits(data, "sf") & (length(variables) == 1)) {
  #   stop('data must be an sf object')
  # }
  #
  if (length(variables) < 1) {
    stop('There should be 1 or more numeric variables')
  }

  data <- data[, variables, drop = FALSE]
  raw_nrow <- nrow(data)
  myNArows <- apply(sf::st_drop_geometry(data), 1, function(x) {
    any(is.na(x))
  })


  data <- stats::na.omit(data)
  data_clust <- data

  if (inherits(data_clust, "sf")) {
    data_clust <-  sf::st_drop_geometry(data_clust)
  }

  my_results <- make_clasification(data_clust,
                                   number_cluster,
                                   fuzzyness = fuzzyness,
                                   distance = distance)


  cluster_na <- data.frame(matrix(NA,
                                  nrow = raw_nrow,
                                  ncol = ncol(my_results$cluster)))
  colnames(cluster_na) <- colnames(my_results$cluster)
  cluster_na[!myNArows, ] <- my_results$cluster
  # Return cluster as character
  cluster_na <- apply(cluster_na, 2, as.character)
  my_results$cluster <- cluster_na
  my_results

}

