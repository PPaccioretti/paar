#' Fuzzy k-means clustering (non-spatial)
#'
#' @description
#' Performs fuzzy k-means clustering on tabular data (non-spatial).
#' This function is a lightweight wrapper around \code{e1071::cmeans},
#' providing a vectorized workflow and clustering quality indices.
#'
#' It is primarily intended as a fallback method when spatial clustering
#' (e.g., \code{kmspc}) cannot be applied, such as when only one variable
#' is available.
#'
#'
#' @inheritParams kmspc
#' @details
#' Missing values are removed prior to clustering. Observations with missing
#' values are reintroduced in the output with \code{NA} cluster assignments.
#'
#' Clustering is performed for each value in \code{number_cluster}, and
#' several indices are returned to assist in selecting the optimal number
#' of clusters:
#' \itemize{
#'   \item Xie-Beni index
#'   \item Partition coefficient
#'   \item Partition entropy
#'   \item Summary index
#' }
#'
#' @return
#' A list with:
#' \describe{
#'   \item{cluster}{\code{data.frame} with cluster assignments for each
#'   evaluated number of clusters}
#'   \item{indices}{\code{data.frame} with clustering validity indices}
#'   \item{summaryResults}{\code{data.frame} with clustering metrics}
#' }
#'
#' @seealso \code{\link{kmspc}}
#' @example inst/examples/fuzzy_k_means.R
#' @export
#'
fuzzy_k_means <- function(
  data,
  variables,
  number_cluster = 3:5,
  fuzzyness = 1.2,
  distance = "euclidean"
) {
  if (missing(variables)) {
    myNumVars <-
      unlist(lapply(sf::st_drop_geometry(data), is.numeric))
    if (sum(myNumVars) == 0) {
      stop('No numeric variables found in data')
    }
    warning("Numeric variables will be used for clustering", call. = FALSE)
    variables <- names(sf::st_drop_geometry(data))[myNumVars]
    message(
      paste("Using variables:", paste(variables, collapse = ", "))
    )
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
    data_clust <- sf::st_drop_geometry(data_clust)
  }

  my_results <- make_clasification(
    data_clust,
    number_cluster,
    fuzzyness = fuzzyness,
    distance = distance
  )

  cluster_na <- data.frame(matrix(
    NA,
    nrow = raw_nrow,
    ncol = ncol(my_results$cluster)
  ))
  colnames(cluster_na) <- colnames(my_results$cluster)
  cluster_na[!myNArows, ] <- my_results$cluster
  # Return cluster as character
  cluster_na <- apply(cluster_na, 2, as.character)
  my_results$cluster <- cluster_na
  my_results
}
