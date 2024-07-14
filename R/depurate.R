#' Remove errors from spatial data
#'
#' @description Data can be filtered by null, edge values, global
#'  outliers and spatial outliers or local defective observations. Default
#'  values are optimized for precision agricultural data.
#'
#' @param x an \code{sf} points object
#' @param y \code{character} with the name of the variable to use for
#'   depuration/filtering process
#' @param toremove \code{character} vector specifying the procedure to
#'   implement for errors removal. Default 'edges', 'outlier', 'inlier'.
#'   See Details.
#' @param crs coordinate reference system: integer with the EPSG code,
#'   or character with proj4string to convert coordinates if \code{x} has
#'   longitude/latitude data
#' @param buffer \code{numeric} distance in meters to be removed. Negative
#'   values are recommended
#' @param ylimitmax \code{numeric} of length 1 indicating the maximum limit
#'   for the \code{y} variable. If \code{NA} \code{Inf} is assumed
#' @param ylimitmin \code{numeric} of length 1 indicating the minimum limit
#'   for the \code{y} variable. If \code{NA} \code{-Inf} is assumed
#' @param sdout \code{numeric} values outside the interval
#'  \eqn{mean ± sdout × sdout} values will be removed
#' @param ldist \code{numeric} lower distance bound to identify neighbors
#' @param udist \code{numeric} upper distance bound to identify neighbors
#' @param criteria \code{character} with "LM" and/or "MP" for methods to
#'    identify spatial outliers
#' @param zero.policy default NULL, use global option value;
#'   if FALSE stop with error for any empty neighbors sets,
#'   if TRUE permit the weights list to be formed with zero-length
#'   weights vectors
#' @param poly_border \code{sf} object with one polygon or NULL. Can be
#' the result of \code{concaveman::concaveman}
#'
#' @details
#' Possible values for \code{toremove} are one or more elements of:
#' \describe{
#'   \item{edges}{All data points for a distance of \code{buffer} m from data
#'   edges are deleted.}
#'   \item{outlier}{Values that are outside the mean±\code{sdout} are removed}
#'   \item{inlier}{Local Moran index of spatial autocorrelation is calculated
#'   for each datum as a tool to identify inliers}
#' }
#'
#' @references
#' Vega, A., Córdoba, M., Castro-Franco, M. et al. Protocol for
#' automating error removal from yield maps. Precision Agric 20, 1030–1044
#' (2019). https://doi.org/10.1007/s11119-018-09632-8
#'
#' @return an object of class \code{paar} with two elements:
#' \describe{
#'  \item{depurated_data}{\code{sf} object with the data after the removal
#'  process}
#'  \item{condition}{\code{character} vector with the condition of each
#'  observation}
#'  }
#' @export
#' @example inst/examples/depurate.R

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
                     criteria = c("LM", "MP"),
                     zero.policy = NULL,
                     poly_border = NULL) {

  toremove <- match.arg(toremove,
                        c('edges', 'outlier', 'inlier'),
                        several.ok = TRUE)

  if (!inherits(x, "sf")) {
    stop('x must be an sf object')
  }
  if (missing(y)) {
    if (ncol(sf::st_drop_geometry(x)) == 1) {
      y <- names(x)[1]
    } else {
      stop('y is missing; must be a valid column name')
    }
  }

  stopifnot('y must be a valid columname' =
              y %in% colnames(x))


  criteria <- match.arg(criteria,
                        c('LM', 'MP'),
                        several.ok = TRUE)


  is_error <- data.frame('idx' = seq_len(nrow(x)),
                         'because' = NA_character_)



  if ('edges' %in% toremove) {
    without_border <- remove_border(
      x = x,
      crs = crs,
      buffer = buffer,
      poly_border = poly_border
    )
    is_error <- is_error_update(is_error,
                                without_border)

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

    x <- without_outlier$depurated_data

  }

  if ('inlier' %in% toremove) {
    without_inlier <- remove_inlier(
      x = x,
      y = y,
      ldist = ldist,
      udist = udist,
      criteria = criteria,
      zero.policy = zero.policy
    )

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
#' @param poly_border \code{sf} object with one polygon or NULL. Can be
#' the result of \code{concaveman::concaveman}
#' @noRd
#' @keywords internal
remove_border <- function(x,
                          crs = NULL,
                          buffer,
                          poly_border = NULL) {
  if (!is.numeric(buffer)) {
    stop("buffer must be numeric", call. = FALSE)
  }

  # Checks if x has longlat crs and change it
  if (is.na(sf::st_crs(x))) {
    warning('Please check results due to crs of object x is NA.')

  }

  if (sf::st_is_longlat(x) & !is.na(sf::st_crs(x))) {
    stopifnot(
      'coordinates provided in crs are longlat degrees' =
        !sf::st_is_longlat(crs),
      'crs must be provided' =
        !is.null(crs)
    )
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

  if (!is.null(poly_border)) {
    if (!inherits(poly_border, "sf")) {
      stop('poly_border must be an sf object', call. = FALSE)
    }

    if (!unique(sf::st_geometry_type(poly_border)) %in% c('POLYGON', 'MULTIPOLYGON')) {
      stop('poly_border must have only POLYGON or MULTIPOLYGON as geometry type',
           call. = FALSE)
    }

    if (nrow(poly_border) != 1) {
      stop(paste(
        'poly_border must have only one POLYGON (1 row)\n',
        'has',
        nrow(poly_border),
        'rows'
      ),
      call. = FALSE)
    }


    mapa_hull <- poly_border

  }


  if (is.null(poly_border)) {
    if (requireNamespace('concaveman', quietly = TRUE)) {
      mapa_hull <- concaveman::concaveman(x)
      message(
        paste0(
          'Concave hull algorithm is computed with\n',
          'concavity = 2 and length_threshold = 0'
        )
      )
    } else {
      message(
        paste0(
          'concaveman package is suggested for this procedure\n',
          'You can install it running install.packages("concaveman")',
          'A convex hull using sf::st_convex_hull will be computed'
        )
      )
      mapa_hull <- sf::st_union(x)
      mapa_hull <- sf::st_convex_hull(mapa_hull)

    }
  }


  mapa_buffer <- sf::st_buffer(mapa_hull, buffer)
  emptyGeometry <- sf::st_is_empty(mapa_buffer)
  if (all(emptyGeometry)) {
    stop(paste0(
      "'buffer' value (",
      buffer,
      ") is higher than all polygons border lengths"
    ),
    call. = FALSE)

  }

  if (any(emptyGeometry)) {
    mapa_buffer <- mapa_buffer[emptyGeometry, ]
    warning(paste0("Some polygons (", sum(emptyGeometry), ") will be ",
                   "removed from the analysis ",
                   "because 'buffer' value (", buffer,
                   ") is higher than polygon length"),
            call. = FALSE)

  }

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
#' @noRd
#' @keywords internal
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
#' @param criteria \code{character} with "LM" and/or "MP" for methods to
#'    identify spatial outliers
#' @param zero.policy default NULL, use global option value;
#'   if FALSE stop with error for any empty neighbors sets,
#'   if TRUE permit the weights list to be formed with zero-length
#'   weights vectors
#' @noRd
#' @keywords internal
remove_inlier <- function(x,
                          y,
                          ldist = 0,
                          udist = 40,
                          criteria = c("LM", "MP"),
                          zero.policy = NULL
                          ) {
  stopifnot('x must be sf class' = inherits(x, 'sf'),
            is.numeric(ldist),
            is.numeric(udist))

  criteria <- match.arg(criteria,
                        c("LM", "MP"),
                        several.ok = TRUE)

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
      lw <- try(spdep::nb2listw(gri, style = 'W', zero.policy = zero.policy),
                silent = TRUE)

      if (all(class(lw) != 'try-error') | i >= maxiter)
        lw_found <- FALSE
      break
    }
  }

  if (!lw_found) {
    stop('Neighbours cannot be identified\n',
         'try modifying the ldist or udist values',
         call. = FALSE)
  }

  condition <- rep(NA_character_, nrow(x))

  # By default non are extracted
  Influ_LM <- rep(FALSE, nrow(x))
  Influ_MP <- rep(FALSE, nrow(x))
  if ("LM" %in% criteria) {
  LM <- spdep::localmoran(x[[y]],
                          lw,
                          # p.adjust.method = 'bonferroni',
                          alternative = 'less',
                          zero.policy = zero.policy)
  # Influence by Local Moran
  # Search column name which start with Pr
  # since spdep >= 1.1.11 change old names
  myColnameProb <- colnames(LM)[grepl("Pr\\(z" , colnames(LM))]

  Influ_LM <- LM[, 'Ii'] < 0 & LM[, myColnameProb] < 0.05
  condition[Influ_LM] <- 'spatial outlier LM'
  }
  # Influence by MoranPlot
  if ("MP" %in% criteria) {
  MP <- spdep::moran.plot(x[[y]],
                          lw,
                          quiet = TRUE,
                          plot = FALSE,
                          zero.policy = zero.policy)

  Influ_MP <- MP$is_inf
  condition[Influ_MP] <- 'spatial outlier MP'
  }
  # x[[condition_name]] <- condition

  mapa_dep <- subset(x, !(Influ_LM | Influ_MP))

  list('depurated_data' = mapa_dep,
       'condition' = condition)
}


#' Auxiliary function for update is_error object in depurate function
#'
#' @param is_error internal object
#' @param remove_result internal results from remove_* functions
#' @noRd
#' @keywords internal
is_error_update <- function(is_error, remove_result) {
  # Keeps NA, others are conditions
  is_error_no_na <- is_error[is.na(is_error$because),]
  idx_remove_result <- which(!is.na(remove_result$condition))

  mycondition <-
    remove_result$condition[!is.na(remove_result$condition)]

  is_error_no_na[idx_remove_result, 'because'] <- mycondition
  # Merge original conditions with new conditions
  is_error_merged <-
    merge(is_error,
          is_error_no_na,
          by = "idx",
          all.x = TRUE)
  # Keeps all conditions of removal
  is_error_merged$because <-
    do.call(pmax, c(is_error_merged[, -1], na.rm = TRUE))

  # Find dots in colnames and remove that columns (they are from merge function)
  is_error_merged[, !agrepl("\\.", colnames(is_error_merged))]

}
