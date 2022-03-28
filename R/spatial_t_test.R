#' Modified t test
#'
#'
#' @description Performs a modified version of the t test to assess the
#' correlation between spatial processes. See SpatialPack::modified.ttest for details.
#'
#' @param data \code{sf} data to extract coordinates or two
#' columns \code{matrix} or \code{data.frame} specifying coordinates.
#'
#' @param variables \code{character} vector with column names to perform ttest
#'
#' @export


spatial_t_test <- function(data, variables) {

  if (!requireNamespace("SpatialPack", quietly = TRUE)) {
    stop(paste0("SpatialPack package is needed to perform spatial_t_test", "\n",
               "Install it with 'install.packages(SpatialPack)'"),
               .call = FALSE)
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
  for (i in 1:n)
  {
    for (j in i:n)
    {
      if (i == j) {
        next

      }

      test <- SpatialPack::modified.ttest(df[, i], df[, j], coords)

      myResults <-
        rbind(myResults,
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
