#' Data check for UV dose calculation
#'
#' UV dose calculation requires date, latitude, and longitude inputs.
#'
#' @param data Input dataframe to be used for UV dose calculation
#' @return A message checking validity of input data
#'
#' @examples
#' mysample <- data.frame(id = c("id000016", "id000021"),
#'       date = as.Date(c("2009-05-15", "2008-08-04")),
#'       easting = c(519000, 365000),
#'       northing = c(176000, 172000))
#' uv_data_check(mysample)
#' @return None, returns a message about input data validity.

#' @export

#check if input data is a dataframe, has date and coordinate columns, has OSGB coordinates
uv_data_check <- function(data) {
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }

  easting_northing <- c("easting", "northing", "east", "north")
  has_osgb <- any(tolower(names(data)) %in% easting_northing)

  date_columns <- sapply(data, lubridate::is.Date)
  has_date_column <- TRUE %in% date_columns

  coordinate_columns <- c("latitude", "longitude", "lat", "lon")
  has_coordinate_columns <- any(tolower(names(data)) %in% coordinate_columns)

    if (has_osgb && !has_coordinate_columns) {
    return("Looks like your data has OSGB coordinates! Use the longitude and latitude functions first to convert 'easting' and 'northing' columns.")
  }

  # Produce error if either date or coordinate columns are not present
  if (!has_date_column || !has_coordinate_columns) {
    stop("Error: The data does not contain a date column and/or coordinate columns.")
  }

  # Return success message
  return("Success! Data contains valid date and coordinate columns.")
}


