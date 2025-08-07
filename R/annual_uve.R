#' Annual erythemal UV
#'
#' This function extracts mean or cumulative erythemal UV dose  at a given date and geographical location (longitude and latitude), over one year prior to the index date.
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @param data data frame containing sample IDs, date, longitude, and latitude
#' @param date A date vector, usually date of assessment or recruitment.
#' @param longitude A numeric vector of longitude values.
#' @param latitude A numeric vector of latitude values.
#' @param temis_path Path to TEMIS UV files downloaded using [temis_uvdec]. Default is current directory.
#' @param region Region of TEMIS data required, options are "europe" (default) or "world".
#' @param type Type of annual dose to be calculated, cumulative (default) or mean.
#' @return A numeric vector of mean or cumulative annual erythemal UV dose estimates.
#' @examples
#' #uses sample TEMIS file
#' mysample <- data.frame(id = c("id000016"),
#'         date = as.Date(c("2010-08-04")),
#'         longitude = c(-2.10),
#'         latitude = c(50.5))
#' uve_example <- system.file("extdata", "uve_example", package="UVdose")
#' annual_uve(mysample, date, longitude, latitude, temis_path=uve_example)
#'
#' @export


annual_uve <- function(data, date, longitude, latitude, temis_path=getwd(), type="cumulative", region="europe"){
  col_date <- substitute(date)
  col_longitude <- substitute(longitude)
  col_latitude <- substitute(latitude)

  val_date <- eval(col_date, data)
  val_longitude <- eval(col_longitude, data)
  val_latitude <- eval(col_latitude, data)

  # Check if date is of class Date
  if (!inherits(val_date, "Date")) {
    stop("Date input must be of class Date.")
  }

  # Check if date is after 2004 from TEMIS data availability
  if (min(val_date) < as.Date("2004-01-19")) {
    stop("Cloud-adjusted UV data is not available prior to 2004-01-19, please remove all earlier index dates from the input.")
  }

  #for non-daily doses it depends on the length of the period required
  if (min(val_date) < as.Date("2004-12-31")) {
    warning("Cloud-adjusted UV data is not available prior to 2004-01-19, please adjust input dates accordingly.")
  }

  # Check if longitude and latitude are numeric
  if (!is.numeric(val_longitude) || !is.numeric(val_latitude)) {
    stop("Longitude and latitude inputs must be numeric.")
  }

  # Check if all inputs are of the same length
  if (!(length(val_date) == length(val_longitude) &&
        length(val_longitude) == length(val_latitude))) {
    stop("All inputs must be of the same length.")
  }

  # Check if any input contains missing data
  if ( anyNA(val_date) || anyNA(val_longitude) || anyNA(val_latitude)) {
    stop("Inputs cannot contain missing values. Remove observations with missing data and try again!")
  }

  # Check if region input is correct
  if (!region %in% c("europe", "world")) {
    stop("Region must be `europe` or `world`.")
  }

  #### CLIMATOLOGY FOR MISSING DATA ####
  #fill in missing data using climatology file
  clim_nc <- ncdf4::nc_open(paste(temis_path,"/", region, "_uvdec_climatology.nc", sep = ""))
  clim_erythema <- ncdf4::ncvar_get(clim_nc, "PRODUCT/uvd_cloudy_mean")
  clim_lat <- ncdf4::ncvar_get(clim_nc, "PRODUCT/latitude")
  clim_lon <- ncdf4::ncvar_get(clim_nc, "PRODUCT/longitude")
  ncdf4::nc_close(clim_nc) # close the connection since we're finished
  # set the dimension names and values to the appropriate latitude and longitude values
  dimnames(clim_erythema) <- list(clim_lon, clim_lat)

  #### RANGE OF LATITUDES/LONGITUDES ####
  # find the range of latitudes and longitudes in the data
  nn <- min(val_longitude) #smallest longitude in the data
  xn <- max(val_longitude) #largest longitude in the data

  # Find the closest value smaller than nn
  index_nn <- findInterval(nn, clim_lon)
  closest_smaller_nn <- clim_lon[index_nn]

  # Find the closest value greater than xn
  index_xn <- findInterval(xn, clim_lon)
  closest_greater_xn <- clim_lon[index_xn + 1]
  message("For longitudes ", nn, " and ", xn, " Nearest grid values are: ", closest_smaller_nn, " and ", closest_greater_xn, "\n")

  #repeat for latitude
  nt <- min(val_latitude)
  xt <- max(val_latitude)
  # Find the closest value smaller than nn
  index_nt <- findInterval(nt, clim_lat)
  closest_smaller_nt <- clim_lat[index_nt]

  # Find the closest value greater than xn
  index_xt <- findInterval(xt, clim_lat)
  closest_greater_xt <- clim_lat[index_xt + 1]
  message("For latitudes ", nt, " and ", xt, " Nearest grid values are: ", closest_smaller_nt, " and ", closest_greater_xt, "\n")

  #### SUBSET CLIMATOLOGY DATA ####
  #select by [lon, lat, day]
  #narrow down the erythema uv array size by limiting to area (lon/lat) of interest
  # longitude [nc 2], latitude [nc 1]

  clim_erythema <- clim_erythema[which(clim_lon == closest_smaller_nn):which(clim_lon == closest_greater_xn),
                                 which(clim_lat == closest_smaller_nt):which(clim_lat == closest_greater_xt), ]

  m_clim_erythema <- reshape2::melt(clim_erythema)
  colnames(m_clim_erythema) <- c("lon", "lat", "day", "uvd")
  clim_erythema <-  reshape2::dcast(m_clim_erythema, lon+lat~day, value.var = "uvd")

  #is this needed?
  clim_erythema <- dplyr::as_tibble(clim_erythema)

  #conflict with 'mutate' from plyr and dplyr
  clim_erythema_leap <- dplyr::mutate(clim_erythema, leap = NA, .after= .data$`59`)
  colnames(clim_erythema_leap) <- c("lon", "lat", 1:366)

  #leap years
  leap <- c(2000, 2004, 2008, 2012, 2016, 2020, 2024)

  #select years of interest
  extract_year_range <- function(date_name){
    year_min <- lubridate::year(min(date_name, na.rm = T)) - 1
    year_max <- lubridate::year(max(date_name, na.rm = T))
    year_range <- year_min:year_max
    return(year_range)
  }

  year <- extract_year_range(val_date)

  # Create a pattern to match the year range
  year_pattern <- paste0("(", paste(year, collapse = "|"), ")")

  #### DAILY UV DATA ####
  #make a list of the nc files downloaded and initiate an empty list
  myfiles <- list.files(path = temis_path, pattern = paste("*uvdec_", region, ".nc", sep=""))

  # Filter files matching the year range
  myfiles <- myfiles[grepl(year_pattern, myfiles)]

  erythema <- list()

  for (i in 1:length(myfiles)){
    # open a conneciton to the ith nc file
    nc_tmp <- ncdf4::nc_open(paste(temis_path, "/", myfiles[i], sep=""))
    # store values from variables and atributes
    nc_erythema <- ncdf4::ncvar_get(nc_tmp, "PRODUCT/uvd_cloudy")[which(clim_lon == closest_smaller_nn):which(clim_lon == closest_greater_xn),
                                                                              which(clim_lat == closest_smaller_nt):which(clim_lat == closest_greater_xt), ]
    nc_lat <- ncdf4::ncvar_get(nc_tmp, "PRODUCT/latitude")[which(clim_lat == closest_smaller_nt):which(clim_lat == closest_greater_xt)]
    nc_lon <- ncdf4::ncvar_get(nc_tmp, "PRODUCT/longitude")[which(clim_lon == closest_smaller_nn):which(clim_lon == closest_greater_xn)]
    # close the connection since we're finished
    ncdf4::nc_close(nc_tmp)
    # set the dimension names and values of your matrix to the appropriate latitude and longitude values
    dimnames(nc_erythema) <- list(lon=nc_lon, lat=nc_lat, day=1:dim(nc_erythema)[3])

    erythema <- append(erythema, list(nc_erythema))
  }

  #name each element of the erythema list after the year
  names(erythema) <- year

  #### COORDINATE MATCHING ####

  #convert the coordinates in mydata data to match the coordinates from the erythema uv data
  #TEMIS coordinates are center points of a range
  #each uv measurement is for a 0.25 latitude x 0.25 longitude square of area
  mylon <- data.frame(nc_lon, lower_lon = c(nc_lon - 0.125), upper_lon = c(nc_lon + 0.125))
  mylat <- data.frame(nc_lat, lower_lat = c(nc_lat - 0.125), upper_lat = c(nc_lat + 0.125))

  data_tmp <- data.frame(date=val_date,longitude=val_longitude, latitude=val_latitude)

  #mydata latitude 57.50000 is exactly at the bound so it gets assigned c(57.375, 57.625)
  #so we use [x >= lower & x < upper] rather than [x >= lower & x <= upper]
  data_tmp$latcoord <- sapply(val_latitude, function(x) {
    out <- mylat$nc_lat[x >= mylat$lower_lat & x < mylat$upper_lat]
    if (length(out) == 0) "NA" else out
  })

  data_tmp$loncoord <- sapply(val_longitude, function(x) {
    out <- mylon$nc_lon[x >= mylon$lower_lon & x < mylon$upper_lon]
    if (length(out) == 0) "NA" else out
  })

  #### RESHAPE THE DATA ####
  #reshape list of 3d arrays into a dataframe then tibble
  merythema <- reshape2::melt(erythema, value.name="uvd")
  colnames(merythema)[which(colnames(merythema) == "L1")] <- "year"

  #is this needed?
  merythema <- dplyr::as_tibble(merythema)

  merythema$date <- as.Date(paste(merythema$year, merythema$day, sep="-"), format="%Y-%j")

  #reduce the size of the UV data, keep only the combinations of lon & lat in the sample data
  data_tmp$coords <- paste(data_tmp$loncoord, data_tmp$latcoord, sep = ",")
  merythema$coords <- paste(merythema$lon, merythema$lat, sep = ",")

  merythema <- dplyr::filter(merythema, .data$coords %in% data_tmp$coords)

  #### FILL IN CLIMATOLOGY ####
  #if uvd in merythema is NA replace with value from climatology file
  clim_erythema <- tidyr::pivot_longer(clim_erythema, cols = !(.data$lon | .data$lat), names_to = "day", values_to = "uvd")
  clim_erythema$day <- as.integer(clim_erythema$day)

  clim_erythema_leap <- tidyr::pivot_longer(clim_erythema_leap, cols = !(.data$lon | .data$lat), names_to = "day", values_to = "uvd")
  clim_erythema_leap$day <- as.integer(clim_erythema_leap$day)

  merythema_clim <- dplyr::left_join(dplyr::filter(merythema, !(year %in% leap)), clim_erythema, by=c("lon", "lat", "day"))
  merythema_clim_leap <- dplyr::left_join(dplyr::filter(merythema, year %in% leap), clim_erythema_leap, by=c("lon", "lat", "day"))

  merythema_clim_all <- dplyr::bind_rows(merythema_clim, merythema_clim_leap) %>%
    dplyr::mutate(uvd.x = dplyr::coalesce(.data$uvd.x, .data$uvd.y))

  #### cumulative or average annual erythemal UV ####
  merythema_clim_all <- dplyr::arrange(merythema_clim_all, .data$coords, .data$date)

  #choose summary function based on input type: cumulative or mean
  summary_fun <- switch(type,
                        "mean" = mean,
                        "cumulative" = sum,
                        stop("Annual UV dose type must be 'cumulative' or 'mean'.")
  )

  #function to calculate annual dose
  uv_year <- function(coords_val, date_val) {
    n <- which(merythema_clim_all$coords == coords_val & merythema_clim_all$date == date_val)
    if (length(n) == 0 || n - 364 < 1) return(NA_real_)
    uv_vals <- merythema_clim_all[(n):(n - 364), "uvd.x", drop = TRUE]
    summary_fun(uv_vals, na.rm = TRUE)
  }

  #apply uv_dose over rows
  ery <- purrr::pmap_dbl(list(data_tmp$coords, data_tmp$date), uv_year)

  return(ery)
}
