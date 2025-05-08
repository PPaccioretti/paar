#' Plot paar data with condition
#'
#' Plots an `paar` data, resulting from `cbind(paar, data)` object colored by a `condition` column.
#' The function supports both `ggplot2` and base R graphics.
#'
#' @param data An `sf` object containing a column named `condition`.
#' @param use_ggplot Logical. If `TRUE` (default), uses `ggplot2` for plotting
#'   if the package is installed. If `FALSE`, or if `ggplot2` is not available,
#'   uses base R plotting.
#'
#' @return If `ggplot2` is used, returns a `ggplot` object invisibly. Otherwise returns `NULL`.
#'
#' @details
#' This function is an S3 method for objects of class `paarCondition`.
#' The input `data` should be an `sf` object with a `condition` column.
#' Typically, this object results from combining filtered result data
#' from the function \link[paar]{depurate},
#' using `cbind(paar_data, data)`, where `paar_data` contains the results
#' returned by `depurate`, and `data` is the original dataset used in the procedure.
#'
#' @export
plot.paarCondition <- function(data, use_ggplot = TRUE, ...) {
  if (!inherits(data, "sf")) {
    stop("The input 'data' must be an sf object.")
  }

  if (!"condition" %in% names(data)) {
    stop("The input data must have a 'condition' column.")
  }

  if (!is.factor(data$condition) && !is.character(data$condition)) {
    warning("'condition' is not a factor or character. Coercing to factor.")
    data$condition <- as.factor(data$condition)
  }

  if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
    condition <- NULL
    p <- ggplot2::ggplot(data) +
      ggplot2::geom_sf(ggplot2::aes(color = condition)) +
      ggplot2::scale_color_discrete(
        labels = function(k) {
          k[is.na(k)] <- "normal"
          k
        },
        na.value = "#44214234"
      ) +
      ggplot2::theme_minimal()
    print(p)
    invisible(p)
  } else {
    condition_levels <- levels(as.factor(data$condition))
    colors <- seq_along(condition_levels)
    plot(data[, "condition"], col = as.numeric(as.factor(data$condition)))
    graphics::legend("topright", legend = condition_levels, fill = colors)
    invisible(NULL)
  }
}
