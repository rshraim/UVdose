
<!-- README.md is generated from README.Rmd. Please edit that file -->

# UVdose

<!-- badges: start -->

[![R-CMD-check](https://github.com/rshraim/UVdose/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rshraim/UVdose/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

UVdose is an R package for the manipulation of UV data from TEMIS. TEMIS
(Tropospheric Emission Monitoring Internet Service) is a web-based
service to browse and download atmospheric satellite data products. More
information about TEMIS is available on their website www.temis.nl.
UVdose was built to facilitate the integration of UV data into health
research. The package functions allow the estimation of daily, seasonal,
and annual doses of erythemal UV and UVB based on geographical
coordinates and dates. Additionally, a function to estimate a cumulative
and weighted UVB dose relevant for vitamin D production in the skin is
available.

## Installation

Install UVdose from [GitHub](https://github.com/) with:

``` r
devtools::install_github("rshraim/UVdose")
```

## Example

Basic usage example:

``` r
library(UVdose)

mysample <- data.frame(id = c("id001"),
         date = as.Date(c("2010-08-04")),
         longitude = c(-2.10),
         latitude = c(50.5))

uvb_example <- system.file("extdata", "uvb_example", package="UVdose")

daily_uvb(mysample, date, longitude, latitude, temis_path = uvb_example)
```
