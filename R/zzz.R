#' @importFrom utils "packageVersion"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(pkgname, " version ", utils::packageVersion(pkgname),
                        paste0("\n", "Thanks for using UVdose! This packages uses data from TEMIS.", "\n",
                               "Please inform TEMIS on publication or reuse of data at <www.temis.nl/contact.php>.", "\n",
                               "Copyright KNMI/ESA <www.temis.nl>"
                      # "Use `citation(\"", pkgname, "\")` to cite in publications. ",
                       ))
}


