#' Compute meridional arc and initial latitude for OSGB conversion
#' @param PHI An integer
#' @param North A numeric vector
#' @noRd


#MARC function: compute meridional arc
Marc <- function(bf0, n, PHI0, PHI){
  M <- bf0 * (((1 + n + ((5 / 4) * (n ^ 2)) + ((5 / 4) * (n ^ 3))) * (PHI - PHI0))
              - (((3 * n) + (3 * (n ^ 2)) + ((21 / 8) * (n ^ 3))) * (sin(PHI - PHI0)) * (cos(PHI + PHI0)))
              + ((((15 / 8) * (n ^ 2)) + ((15 / 8) * (n ^ 3))) * (sin(2 * (PHI - PHI0))) * (cos(2 * (PHI + PHI0))))
              - (((35 / 24) * (n ^ 3)) * (sin(3 * (PHI - PHI0))) * (cos(3 * (PHI + PHI0)))))
  return(M)
}

#InitialLat function: Compute initial value for Latitude (PHI) IN RADIANS.
InitialLat <- function(North, n0, af0, PHI0, n, bf0){
  #First PHI value (PHI1)
  PHI1 <- ((North - n0) / af0) + PHI0
  #Calculate M
  M <- Marc(bf0, n, PHI0, PHI1)
  #Calculate new PHI value (PHI2)
  PHI2 <- ((North - n0 - M) / af0) + PHI1
  #Iterate to get final value for InitialLat
  while (abs(North - n0 - M) > 1E-05){
    PHI2 <- ((North - n0 - M) / af0) + PHI1
    M <- Marc(bf0, n, PHI0, PHI2)
    PHI1 <- PHI2
  }
  return(PHI2)
}

