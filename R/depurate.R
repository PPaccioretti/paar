#' Remove errors from spatial data
#'
#' @description Data can be filtered by null, edge values, global
#'  outliers and spatial outliers or local defective observations. Default
#'  values are optimized for precision agricultural data.
#'
#' @inheritParams remove_border
#' @inheritParams remove_outlier
#' @inheritParams remove_inlier
#' @param toremove \code{character} vector specifying the procedure to
#'   implement for errors removal. Default 'edges', 'outlier', 'inlier'.
#'   See Details.
#'
#' @details Possible values for \code{toremove} are one or more elements of:
#' \describe{
#'   \item{edges}{All data points for a distance of \code{buffer} m from data
#'   edges are deleted}
#'   \item{outlier}{Values that are outside the mean±\code{sdout} are removed}
#'   \item{inlier}{Local Moran index of spatial autocorrelation is calculated
#'   for each datum as a tool to identify inliers}
#' }
#'
#' @references Vega, A., Córdoba, M., Castro-Franco, M. et al. Protocol for
#'   automating error removal from yield maps. Precision Agric 20, 1030–1044
#'   (2019). https://doi.org/10.1007/s11119-018-09632-8
#' @return an object of class \code{paar}
#' @export
#'
#' @examples
depurate <- function(x,
                     y,
                     toremove = c('edges', 'outlier', 'inlier'),
                     crs = NULL,
                     buffer = -10,
                     ylimitmax = NA,
                     ylimitmin = 0,
                     sdout = 3,
                     ldist = 0,
                     udist = 40,
                     zero.policy = NULL) {

  toremove <- match.arg(toremove,
                        c('edges', 'outlier', 'inlier'),
                        several.ok = TRUE)

  stopifnot('y must be a vald columname' =
              y %in% colnames(x))

  is_error <- data.frame('idx' = seq_len(nrow(x)),
                         'because' = NA_character_)



  if ('edges' %in% toremove) {
    without_border <- remove_border(
      x = x,
      crs = crs,
      buffer = buffer
    )
    is_error <- is_error_update(is_error,
                                without_border)
    # is_error[which(!is.na(without_border$condition)), 'because'] <- 'edges'

    x <- without_border$depurated_data
  }

  if ('outlier' %in% toremove) {
    without_outlier <- remove_outlier(
      x = x,
      y = y,
      ylimitmax = ylimitmax,
      ylimitmin = ylimitmin,
      sdout = sdout
    )
    is_error <- is_error_update(is_error,
                                without_outlier)
    # is_error[is.na(is_error$because), ][which(!is.na(without_outlier$condition)),
    #                                     'because'] <-
    #   without_outlier$condition[!is.na(without_outlier$condition)]

    x <- without_outlier$depurated_data

  }

  if ('inlier' %in% toremove) {
    without_inlier <- remove_inlier(
      x = x,
      y = y,
      ldist = ldist,
      udist = udist,
      zero.policy = zero.policy
    )

    # is_error[is.na(is_error$because), ][which(!is.na(without_inlier$condition)),
    #                                     'because'] <-
    #   without_inlier$condition[!is.na(without_inlier$condition)]
    is_error <- is_error_update(is_error,
                                without_inlier)
    x <- without_inlier$depurated_data
  }

  results <- list('depurated_data' = x,
                  'condition' = is_error$because)

  results <- structure(results,
                       class = c('list', 'paar'))
  results

}




#' Remove borders
#'
#' @details Removes all points from \code{x} that are \code{buffer} meters from
#'   boundary.
#'
#' @inheritParams remove_inlier
#' @param crs coordinate reference system: integer with the EPSG code,
#'   or character with proj4string to convert coordinates if \code{x} has
#'   longitude/latitude data
#' @param buffer \code{numeric} distance in meters to be removed. Negative
#'   values are recommended
#'
remove_border <- function(x,
                          crs = NULL,
                          buffer) {
  stopifnot(is.numeric(buffer))

  # Checks if x has longlat crs and change it
  if (is.na(sf::st_crs(x))) {
    warning('Pleas check results due to crs of object x is NA.')

  }

  if (sf::st_is_longlat(x) & !is.na(sf::st_crs(x))) {
    stopifnot('coordinates provided in crs are longlat degrees' =
                !sf::st_is_longlat(crs),
              'crs must be provided' =
                !is.null(crs))
    x <- sf::st_transform(x, crs = crs)

  }

  # Checks if units package is installed to make buffer in meters

    if (!is.na(sf::st_crs(x))) {
      if (requireNamespace('units', quietly = TRUE)) {
        buffer <- units::as_units(buffer, 'm')
      } else {
        message(
          paste0(
            'units package is suggested for this procedure\n',
            'You can install it running install.packages("units")'
          )
        )
      }
    }

  if (as.numeric(buffer) > 0) {
    message('Negative buffer value is recommended')
  }

  # mapa_hull <- sf::st_union(x)
  # mapa_hull <- sf::st_convex_hull(mapa_hull)
  mapa_hull <- concaveman::concaveman(x)
  mapa_buffer <- sf::st_buffer(mapa_hull, buffer)
  is_inside <- sf::st_intersects(x, mapa_buffer)
  is_inside <- do.call(c, lapply(is_inside, length))
  is_inside <- as.logical(is_inside)

  mapa_dep <- x[is_inside, ]

  condition <- rep(NA_character_, nrow(x))
  condition[!is_inside] <- 'border'
  # x[[condition_name]] <- condition

  list('depurated_data' = mapa_dep,
       'condition' = condition)

}

#' @title Removes outliers
#'
#' @inheritParams remove_inlier
#' @param ylimitmax \code{numeric} of length 1 indicating the maximum limit
#'   for the \code{y} variable. If \code{NA} \code{Inf} is assumed
#' @param ylimitmin \code{numeric} of length 1 indicating the minimum limit
#'   for the \code{y} variable. If \code{NA} \code{-Inf} is assumed
#' @param sdout \code{numeric} values outside the interval
#'  \eqn{mean ± sdout × sdout} values will be removed
#'
remove_outlier <- function(x,
                           y,
                           ylimitmax = NA,
                           ylimitmin = 0,
                           sdout = 3) {
  stopifnot(
    'ylimitmax must be numeric or NA' =
      is.numeric(ylimitmax) |
      is.na(ylimitmax) | is.null(ylimitmax),
    'ylimitmin must be numeric or NA' =
      is.numeric(ylimitmin) |
      is.na(ylimitmin) | is.null(ylimitmin),
    'y must be a vald columname' =
      y %in% colnames(x)
  )

  if (is.na(ylimitmax) | is.null(ylimitmax)) {
    ylimitmax <- Inf
  }

  if (is.na(ylimitmin) | is.null(ylimitmin)) {
    ylimitmax <- -Inf
  }


  global_min <- x[[y]] <= ylimitmin
  global_max <- x[[y]] >= ylimitmax
  mapa_dep <- subset(x, !global_min & !global_max)

  mean_mapa <- mean(mapa_dep[[y]], na.rm = TRUE)
  stdev_mapa <- sqrt(stats::var(mapa_dep[[y]], na.rm = TRUE))

  LI <- mean_mapa - sdout * stdev_mapa
  LS <- mean_mapa + sdout * stdev_mapa

  outlier_dep <- mapa_dep[[y]] <= LI | mapa_dep[[y]] >= LS
  outlier <- x[[y]] <= LI | x[[y]] >= LS

  mapa_dep <- subset(mapa_dep, !outlier_dep)

  condition <- rep(NA_character_, nrow(x))
  condition[outlier] <- 'outlier'
  condition[global_min] <- 'global min'
  condition[global_max] <- 'global max'

  # x[[condition_name]] <- condition

  list('depurated_data' = mapa_dep,
       'condition' = condition)

}

#' @title Remove spatial outliers
#' @description Removes spatial outliers using Local Moran's I statistic
#'   and moran scatterplot.
#'
#' @param x an \code{sf} points object
#' @param y \code{character} with the name of the variable to use for
#'    depuration process
#' @param ldist \code{numeric} lower distance bound to identify neighbors
#' @param udist \code{numeric} upper distance bound to identify neighbors
#' @param zero.policy default NULL, use global option value;
#'   if FALSE stop with error for any empty neighbors sets,
#'   if TRUE permit the weights list to be formed with zero-length
#'   weights vectors
#'
remove_inlier <- function(x,
                          y,
                          ldist = 0,
                          udist = 40,
                          zero.policy = NULL) {
  stopifnot('x must be sf class' = inherits(x, 'sf'),
            is.numeric(ldist),
            is.numeric(udist))

  maxiter <- 10

  gri <- spdep::dnearneigh(x, ldist, udist)
  lw <-
    try(spdep::nb2listw(gri, style = 'W', zero.policy = zero.policy))

  umin <- udist
  i <- 1
  lw_found <- TRUE

  if (inherits(lw, 'try-error')) {
    repeat {
      lw_found <- FALSE
      i <- i + 1
      umin <- umin + 10
      gri <- spdep::dnearneigh(x, ldist, umin)
      lw <- try(spdep::nb2listw(gri, style = 'W'),
                silent = TRUE)

      if (all(class(lw) != 'try-error') | i >= maxiter)
        lw_found <- FALSE
      break
    }
  }

  if (!lw_found) {
    stop('Neighbours cannot be identified\n',
         'try modifying the ldist or udist values')
  }

  LM <- spdep::localmoran(x[[y]],
                          lw,
                          p.adjust.method = 'bonferroni',
                          alternative = 'less')

  MP <- spdep::moran.plot(x[[y]],
                          lw,
                          quiet = TRUE,
                          plot = FALSE,
                          zero.policy = zero.policy)

  # Influence by Local Moran
  Influ_LM <- LM[, 'Ii'] < 0 & LM[, 'Pr(z < 0)'] < 0.05
  # Influence by MoranPlot
  Influ_MP <- MP$is_inf

  condition <- rep(NA_character_, nrow(x))
  condition[Influ_LM] <- 'spatial outlier LM'
  condition[Influ_MP] <- 'spatial outlier MP'

  # x[[condition_name]] <- condition

  mapa_dep <- subset(x, !(Influ_LM | Influ_MP))

  list('depurated_data' = mapa_dep,
       'condition' = condition)
}


#' Auxiliary function for update is_error object in depurate function
#'
#' @param is_error internal object
#' @param remove_result internal results from remove_* functions
#'
is_error_update <- function(is_error, remove_result) {
  # Keeps NA, others are conditions
  is_error_no_na <- is_error[is.na(is_error$because), ]
  idx_remove_result <- which(!is.na(remove_result$condition))

  mycondition <-
    remove_result$condition[!is.na(remove_result$condition)]

  is_error_no_na[idx_remove_result, 'because'] <- mycondition
  # Merge original conditions with new conditions
  is_error_merged <-
    merge(is_error,
          is_error_no_na,
          by = "idx",
          all.x = TRUE
    )
  # Keeps all conditions of removal
  is_error_merged$because <-
    do.call(pmax, c(is_error_merged[ , -1], na.rm = TRUE))

  # Find dots in colnames and remove that columns (they are from merge function)
  is_error_merged[,!agrepl("\\.", colnames(is_error_merged))]

}

