#' Download TEMIS file
#'
#' These functions download UV files from TEMIS for a given range of years present. `temis_uvdvc` returns vitamin D UV data and `temis_uvdec` returns erythemal UV data. `temis_clim` returns only a climatology file.
#' For various technical reasons, some days are missing UV observations. TEMIS provides a climatology file, which is a UV file with values averaged across 2004-2020 and is used in the other functions to fill in these missing observations.
#' In the first two functions, the climatology file is downloaded by default and used downstream to fill in missing observations. It contains UV (erythemal or UVB) values for each day averaged over 17 years: 2004 - 2020 (leap day 29 Feb. is skipped).
#' Yearly UV files are downloaded from TEMIS for Europe by default.
#' The UV doses, in kJ/m2, are calculated based on cloud-adjusted data. See maps below for area coverage for each of the "europe" and "world" regions (as of 21 March 2025).
#' See <https://www.temis.nl/uvradiation/product/uvncinfo.html> for more info. Cloud-adjusted TEMIS data is available from 19-Jan-2004 onwards.
#' @importFrom utils "download.file"
#' @param years The range of years for which UV files will be downloaded, either an integer range or a date vector such as a date column in a dataframe. The data file for an additional year to the provided range is downloaded for non-daily UV dose calculations to account for earlier dates. For example, for the CW-D-UVB dose in [cw_uvb] a sample dated 01-02-2007 requires UV data up to 18-09-2006.
#' @param path Directory where files will be downloaded to ("path/to/dir").
#' @param climatology If TRUE (default) the climatology file will be downloaded. If FALSE, only year files will be downloaded. The same climatology file is used regardless of the specified year range.
#' @param region As illustrated in the maps below, if "europe" (default), files covering Europe region will be downloaded, if "world", world files will be downloaded.
#' \if{html}{Europe region coverage:\out{<div style="text-align: left">}\figure{March2025_vitD_cloudy_europe.png}{options: style="width:500px;max-width:50\%;"}\out{</div><p>}}
#' \if{latex}{Europe region coverage:\out{\begin{center}}\figure{March2025_vitD_cloudy_europe.png}\out{\end{center}}}
#' \if{html}{World region coverage:\out{<div style="text-align: left">}\figure{March2025_vitD_cloudy_world.png}{options: style="width:500px;max-width:50\%;"}\out{</div>}}
#' \if{latex}{Europe region coverage:\out{\begin{center}}\figure{March2025_vitD_cloudy_world.png}\out{\end{center}}}
#' @param uv_type For temis_clim, "uve" or "uvb" for erythemal UV or vitamin D UVB, respectively.
#' @return Files downloaded to specified directory.


#' @export
# vitamin D UV
temis_uvdvc <-  function(years, path, climatology = TRUE, region = "europe"){
 if(region == "europe"){
   region_nc <- "_europe.nc"
 } else if (region == "world") {
   region_nc <- "_world.nc"
 } else {
   stop("Error in region input. Must be europe or world.")
 }

   if(climatology == TRUE){
    message("Downloading TEMIS Vitamin D UV climatology file")
    download.file(paste("https://www.temis.nl/uvradiation/v2.0/nc/clim/uvdvcclim", region_nc, sep=""),
                  destfile = paste(path, "/",region, "_uvdvc_climatology.nc", sep=""))
  }

  if (is.integer(years)) {
    year_range <- c(min(years)-1, years)
    for (y in year_range){
      download.file(paste("https://www.temis.nl/uvradiation/v2.0/nc/",y,"/uvdvc",y,region_nc, sep = ""),
                    destfile = paste(path, "/", y, "_uvdvc", region_nc, sep = ""))
    }
  } else if (lubridate::is.Date(years)) {
    year_min <- lubridate::year(min(years, na.rm = T)) - 1
    year_max <- lubridate::year(max(years, na.rm = T))
    year_range <- year_min:year_max

    message(paste("Downloading TEMIS annual vitamin D UV files", paste(year_range, collapse = ","), sep=" "))

    for (y in year_range){
      download.file(paste("https://www.temis.nl/uvradiation/v2.0/nc/",y,"/uvdvc",y,region_nc, sep = ""),
                    destfile = paste(path,"/",y, "_uvdvc",region_nc, sep = ""))
    }
  } else {
    stop("Input must be either an integer year range (e.g. 2010:2012) or a date vector.")
  }
  message("Heads up! Do not rename files, downstream UV functions search by file name.")
}


#' @export
#' @rdname temis_uvdvc
# erythemal UV
temis_uvdec <-  function(years, path, climatology = TRUE, region = "europe"){
  if(region == "europe"){
    region_nc <- "_europe.nc"
  } else if (region == "world") {
    region_nc <- "_world.nc"
  } else {
    stop("Error in region input. Must be europe or world.")
  }

  if(climatology == TRUE){
    message("Downloading TEMIS erythemal UV climatology file")
    download.file(paste("https://www.temis.nl/uvradiation/v2.0/nc/clim/uvdecclim", region_nc, sep=""),
                  destfile = paste(path, "/",region, "_uvdec_climatology.nc", sep=""))
  }

  if (is.integer(years)) {
    year_range <- c(min(years)-1, years)
    for (y in year_range){
      download.file(paste("https://www.temis.nl/uvradiation/v2.0/nc/",y,"/uvdec",y,region_nc, sep = ""),
                    destfile = paste(path, "/",y, "_uvdec", region_nc, sep = ""))
    }
  } else if (lubridate::is.Date(years)) {
    year_min <- lubridate::year(min(years, na.rm = T)) - 1
    year_max <- lubridate::year(max(years, na.rm = T))
    year_range <- year_min:year_max

    message(paste("Downloading TEMIS annual erythemal UV files", paste(year_range, collapse = ","), sep=" "))

    for (y in year_range){
      download.file(paste("https://www.temis.nl/uvradiation/v2.0/nc/",y,"/uvdec",y,region_nc, sep = ""),
                    destfile = paste(path,"/",y, "_uvdec", region_nc, sep = ""))
    }
  } else {
    stop("Input must be either an integer year range (e.g. 2010:2012) or a date vector.")
  }
  message("Heads up! Do not rename files, downstream UV functions search by file name.")
}


#' @export
#' @rdname temis_uvdvc
#climatology file only
temis_clim <-  function(path, uv_type, region = "europe"){
  if(region == "europe"){
    region_nc <- "_europe.nc"
  } else if (region == "world") {
    region_nc <- "_world.nc"
  } else {
    stop("Error in region input. Must be europe or world.")
  }

  if (uv_type == "uve"){
    message("Downloading TEMIS erythemal UV climatology file")
    download.file(paste("https://www.temis.nl/uvradiation/v2.0/nc/clim/uvdecclim", region_nc, sep=""),
                  destfile = paste(path, "/",region, "_uvdec_climatology.nc", sep=""))
  } else if (uv_type == "uvb") {
    message("Downloading TEMIS Vitamin D UV climatology file")
    download.file(paste("https://www.temis.nl/uvradiation/v2.0/nc/clim/uvdvcclim", region_nc, sep=""),
                  destfile = paste(path, "/",region, "_uvdvc_climatology.nc", sep=""))
  } else {
    stop("Type input must be either uve (erythemal UV) or uvb (vitamin D UVB).")
  }
  message("Heads up! Do not rename files, downstream UV functions search by file name.")
}




