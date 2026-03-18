#' Modified t test
#'
#'
#' @description
#' Performs a modified t-test to assess the correlation between variables
#' while accounting for spatial autocorrelation. This implementation wraps
#' \code{SpatialPack::modified.ttest}.
#'
#' @details
#' The function computes pairwise correlations between the specified variables
#' and adjusts the significance test to account for spatial dependence using
#' coordinates. If \code{data} is an \code{sf} object, coordinates are extracted
#' automatically. Otherwise, coordinates must be provided as an object with two
#' columns.
#'
#' @param data An \code{sf} object containing geometry and variables, or a
#'   \code{matrix}/\code{data.frame} with two columns representing
#'   spatial coordinates (e.g., X and Y).
#'
#' @param variables A \code{character} vector with the names of the variables
#'   to be tested. If \code{data} is not an \code{sf} object, this should be a
#'   matrix or data.frame of variables to test.
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{Var1}{Name of the first variable}
#'   \item{Var2}{Name of the second variable}
#'   \item{corr}{Estimated correlation coefficient}
#'   \item{p.value}{P-value adjusted for spatial autocorrelation}
#' }
#'
#' @seealso \code{\link[SpatialPack]{modified.ttest}}
#'
#' @example inst/examples/spatial_t_test.R
#'
#' @export

spatial_t_test <- function(data, variables) {
  if (!requireNamespace("SpatialPack", quietly = TRUE)) {
    stop(
      paste0(
        "SpatialPack package is needed to perform spatial_t_test",
        "\n",
        "Install it with 'install.packages(\'SpatialPack\')'"
      ),
      .call = FALSE
    )
  }

  #MethodTitle= Test t modificado
  #MethodNickName= CorrEsp
  if (inherits(data, "sf")) {
    stopifnot(is.character(variables))
    data <- stats::na.omit(data[, variables])
    coords <- sf::st_coordinates(data)
    variables <- sf::st_drop_geometry(data[, variables])
  } else {
    stopifnot(ncol(data) == 2)
    mydata <- as.data.frame(stats::na.omit(cbind(data, variables)))
    coords <- as.matrix(mydata[, 1:2])
    variables <- mydata[, -c(1:2)]
  }

  n <- ncol(variables)
  df <- variables
  myResults <- data.frame()
  for (i in 1:n) {
    for (j in i:n) {
      if (i == j) {
        next
      }

      test <- SpatialPack::modified.ttest(df[, i], df[, j], coords)

      myResults <-
        rbind(
          myResults,
          data.frame(
            "Var1" = names(df)[i],
            "Var2" = names(df)[j],
            "corr" = test$corr,
            p.value = test$p.value
          )
        )
    }
  }

  myResults
}
