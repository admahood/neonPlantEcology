#' Polygons of National Ecological Observatory Network Core and Relocatable Terrestrial Sites
#'
#' Note: Some sites have more than one polygon. There are 59 polygons and 47 total sites.
#'
#' @format ## 'site_polygons'
#' A simple feature collection with 59 features and 8 fields
#' \describe{
#'  \item{domainNumb}{Domain Number}
#'  \item{domainName}{Domain Name}
#'  \item{siteType}{Site type. Core or Relocatable}
#'  \item{siteName}{Site name}
#'  \item{siteID}{Four letter site ID. Used in npe_download()}
#'  \item{siteHost}{Organization hosting the site}
#'  \item{areaKm2}{Area of the site in square kilometers}
#'  \item{acres}{Area of the site in acres}
#'  \item{geometry}{list column containing geometry information for each polygon}
#' }
#' @source <https://www.neonscience.org>
"site_polygons"

#' National Ecological Observatory Network Core and Relocatable Terrestrial Sites
#'
#' Note: Some sites have more than one polygon. There are 59 polygons and 47 total sites.
#'
#' @format ## 'sites'
#' data frame with 59 features and 7 fields
#' \describe{
#'  \item{domainNumb}{Domain Number}
#'  \item{domainName}{Domain Name}
#'  \item{siteType}{Site type. Core or Relocatable}
#'  \item{siteName}{Site name}
#'  \item{siteID}{Four letter site ID. Used in npe_download()}
#'  \item{siteHost}{Organization hosting the site}
#'  \item{areaKm2}{Area of the site in square kilometers}
#'  \item{acres}{Area of the site in acres}
#' }
#' @source <https://www.neonscience.org>
"sites"


#' Plant Presence and Percent Cover Data for Domain 14
#'
#' This includes Jornada Experimental Range and Santa Rita Experimental Range
#'
#' @format ## 'D14'
#' A list with 8 items, 2 of which are used by neonPlantEcology
#' @source <https://doi.org/10.48443/9579-a253>
#' @source <https://data.neonscience.org/data-products/DP1.10058.001/RELEASE-2023>
"D14"

#' Plot centroids for the entire NEON network
#'
#' @format ## 'plot_centroids'
#' A simple feature collection with 3842 features and 36 fields
#' @source <https://www.neonscience.org>
"plot_centroids"
