#' Barley grain yield
#'
#' A dataset containing Barley grain yield using calibrated commercial
#' yield monitors mounted on combines equipped with DGPS.
#'
#' @format A data frame with 7395 rows and 3 variables:
#' \describe{
#'   \item{X}{X coordinate, in meters}
#'   \item{Y}{Y coordinate, in meters}
#'   \item{Yield}{grain yield, in ton per hectare}
#' }
#' @details
#' Coordinate reference system is "WGS 84 / UTM zone 20S", epsg:32720
"barley"


#' Database from a production field under continuous agriculture
#'
#' A database from a wheat (Triticum aestivum L.) production field (60 ha)
#' under
#' continuous agriculture, located in south-eastern Pampas, Argentina.
#'
#' @format A data frame with 5982 rows and 7 variables:
#' \describe{
#'   \item{x}{X coordinate, in meters}
#'   \item{y}{Y coordinate, in meters}
#'   \item{CE30}{apparent electrical conductivity taken at 0–30 cm}
#'   \item{CE90}{apparent electrical conductivity taken at 0–90 cm}
#'   \item{Elev}{elevation, in meters}
#'   \item{Pe}{soil depth, in centimeters}
#'   \item{Tg}{wheat grain yield}
#' }
#' @details
#' Coordinate reference system is "WGS 84 / UTM zone 20S", epsg:32720
#' Wheat grain yield was recorded in 2009 using calibrated commercial yield
#' monitors mounted on combines equipped with DGPS. Soil ECa measurements were taken
#' using Veris 3100 (VERIS technologies enr., Salina, KS, USA). Soil depth was
#' measured using a hydraulic penetrometer on a 30 × 30 m regular
#' grid (Peralta et al., 2015). Re-gridding was performed to obtain values of
#' all variables at each intersection point of a 10 × 10 m grid.
#'
#' @references
#' N.R. Peralta, J.L. Costa, M. Balzarini, M. Castro Franco, M. Córdoba,
#' D. Bullock
#' Delineation of management zones to improve nitrogen management of wheat
#' Comput. Electron. Agric., 110 (2015), pp. 103-113,
#' 10.1016/j.compag.2014.10.017
#'
#' Paccioretti, P., Córdoba, M., & Balzarini, M. (2020). FastMapping:
#' Software to create field maps and identify management zones in precision
#' agriculture. Computers and Electronics in Agriculture, 175, 105556.
"wheat"
