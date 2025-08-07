#' Un-project Transverse Mercator Eastings and Northings back to latitude and longitude.
#'
#' These functions convert OSGB (Ordnance Survey of Great Britain) coordinates, i.e. Transverse Mercator easting and northing coordinates,
#' for example, as provided by the UK Biobank.
#' `latitude` returns latitude values and `longitude` returns longitude values.
#'
#' @param data a data frame containing OSGB coordinates
#' @param easting a numeric vector of easting coordinates
#' @param northing a numeric vector of northing coordinates
#' @examples
#' osgb <- data.frame(east = c(393000, 461000, 438000), north = c(287000, 223000, 565000))
#' latitude(osgb, east, north)
#' longitude(osgb, east, north)
#' @name OSGB
NULL


#' @rdname OSGB
#' @export
latitude <- function(data, easting, northing){

  ##CONSTANTS
  #ellipsoid axis dimensions (a & b) in meters; _
  a <- 6377563.3960 #a for OSGB36
  b <- 6356256.9090 #b for OSGB36
  #eastings (e0) and northings (n0) of false origin in meters; _
  e0 <- 400000.000
  n0 <- -100000.000
  #central meridian scale factor (f0) and _
  f0 <- 0.999601271700
  #latitude (PHI0) and longitude (LAM0) of false origin in decimal degrees.
  PHI0 <- 49 #phi 0 N
  LAM0 <- -2 #lambda 0 W
  #Convert angle measures to radians
  RadPHI0 <- PHI0 * (pi / 180)
  RadLAM0 <- LAM0 * (pi / 180)
  #Compute af0, bf0, e squared (e2), n and Et
  af0 <- a * f0
  bf0 <- b * f0
  e2 <- ((af0 ^ 2) - (bf0 ^ 2)) / (af0 ^ 2)
  n <- (af0 - bf0) / (af0 + bf0)


  col_east <- substitute(easting)
  col_north <- substitute(northing)

  East <- eval(col_east, data)
  North <- eval(col_north, data)

  # Check if longitude and latitude are numeric
  if (!is.numeric(East) || !is.numeric(North)) {
    stop("Longitude and latitude inputs must be numeric.")
  }

  # Check if all inputs are of the same length
  if (length(East) != length(North)) {
    stop("All inputs must be of the same length.")
  }

  # Check if any input contains missing data
  if (anyNA(East) || anyNA(North)) {
    stop("Inputs cannot contain missing values. Remove observations with missing data and try again!")
  }

  # Initialize an empty vector
  latitude <- numeric(length(East))

  for (i in seq_along(East)){
  Et <- East[i] - e0
  #Compute initial value for latitude (PHI) in radians
  PHId <- InitialLat(North[i], n0, af0, RadPHI0, n, bf0)
  #Compute nu, rho and eta2 using value for PHId
  nu <- af0 / (sqrt(1 - (e2 * ((sin(PHId)) ^ 2))))
  rho <- (nu * (1 - e2)) / (1 - (e2 * (sin(PHId)) ^ 2))
  eta2 <- (nu / rho) - 1
  ##LATITUDE
  #Compute Latitude
  VII <- (tan(PHId)) / (2 * rho * nu)
  VIII <- ((tan(PHId)) / (24 * rho * (nu ^ 3))) * (5 + (3 * ((tan(PHId)) ^ 2)) + eta2 - (9 * eta2 * ((tan(PHId)) ^ 2)))
  IX <- ((tan(PHId)) / (720 * rho * (nu ^ 5))) * (61 + (90 * ((tan(PHId)) ^ 2)) + (45 * ((tan(PHId)) ^ 4)))
  latitude[i] <- (180 / pi) * (PHId - ((Et ^ 2) * VII) + ((Et ^ 4) * VIII) - ((Et ^ 6) * IX))

  }
  return(latitude)
}

#' @rdname OSGB
#' @export
longitude <- function(data, easting, northing){


  ##CONSTANTS
  #ellipsoid axis dimensions (a & b) in meters; _
  a <- 6377563.3960 #a for OSGB36
  b <- 6356256.9090 #b for OSGB36
  #eastings (e0) and northings (n0) of false origin in meters; _
  e0 <- 400000.000
  n0 <- -100000.000
  #central meridian scale factor (f0) and _
  f0 <- 0.999601271700
  #latitude (PHI0) and longitude (LAM0) of false origin in decimal degrees.
  PHI0 <- 49 #phi 0 N
  LAM0 <- -2 #lambda 0 W
  #Convert angle measures to radians
  RadPHI0 <- PHI0 * (pi / 180)
  RadLAM0 <- LAM0 * (pi / 180)
  #Compute af0, bf0, e squared (e2), n and Et
  af0 <- a * f0
  bf0 <- b * f0
  e2 <- ((af0 ^ 2) - (bf0 ^ 2)) / (af0 ^ 2)
  n <- (af0 - bf0) / (af0 + bf0)

  col_east <- substitute(easting)
  col_north <- substitute(northing)

  East <- eval(col_east, data)
  North <- eval(col_north, data)

  # Check if longitude and latitude are numeric
  if (!is.numeric(East) || !is.numeric(North)) {
    stop("Longitude and latitude inputs must be numeric.")
  }

  # Check if all inputs are of the same length
  if (length(East) != length(North)) {
    stop("All inputs must be of the same length.")
  }

  # Check if any input contains missing data
  if (anyNA(East) || anyNA(North)) {
    stop("Inputs cannot contain missing values. Remove observations with missing data and try again!")
  }

  # Initialize an empty vector
  longitude <- numeric(length(East))
  for (i in seq_along(East)){
  Et <- East[i] - e0
  #Compute initial value for latitude (PHI) in radians
  PHId <- InitialLat(North[i], n0, af0, RadPHI0, n, bf0)
  #Compute nu, rho and eta2 using value for PHId
  nu <- af0 / (sqrt(1 - (e2 * ((sin(PHId)) ^ 2))))
  rho <- (nu * (1 - e2)) / (1 - (e2 * (sin(PHId)) ^ 2))
  eta2 <- (nu / rho) - 1
  ##LONGITUDE
  #Compute Longitude
  X <- ((cos(PHId)) ^ -1) / nu
  XI <- (((cos(PHId)) ^ -1) / (6 * (nu ^ 3))) * ((nu / rho) + (2 * ((tan(PHId)) ^ 2)))
  XII <- (((cos(PHId)) ^ -1) / (120 * (nu ^ 5))) * (5 + (28 * ((tan(PHId)) ^ 2)) + (24 * ((tan(PHId)) ^ 4)))
  XIIA <- (((cos(PHId)) ^ -1) / (5040 * (nu ^ 7))) * (61 + (662 * ((tan(PHId)) ^ 2)) + (1320 * ((tan(PHId)) ^ 4)) + (720 * ((tan(PHId)) ^ 6)))
  longitude[i] <- (180 / pi) * (RadLAM0 + (Et * X) - ((Et ^ 3) * XI) + ((Et ^ 5) * XII) - ((Et ^ 7) * XIIA))

  }
  return(longitude)
}


