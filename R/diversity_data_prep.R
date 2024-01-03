
#' Data downloader
#'
#' A wrapper function to download data from the NEON API. Some commonly used
#' products are provided as plain language options, otherwise the user
#' can enter the product ID number (dpID). Downloads Plant Presence and Percent
#' Cover by default (DP1.10058.001).
#'
#' @param sites a vector of NEON site abbreviations. Defaults to "JORN"
#' @param product a plain language vector of the data product to be downloaded.
#' Can be "plant_diversity", "litterfall", "woody_veg_structure",
#' "belowground_biomass", "herbaceous_clip", "coarse_downed_wood",
#' or "soil_microbe_biomass"
#' @param dpID if you need a data product not given as one of the product
#' options, set the data product ID here (e.g. "DP1.10023.001").
#' @keywords download neon diversity
#'
#' @examples
#' \dontrun{diversity_object <- npe_download(sites = "JORN")}
#' @returns a list
#' @export
npe_download <- function(sites = "JORN",
                         dpID = NA,
                         product = "plant_diversity"){
  requireNamespace("neonUtilities")

  if(!is.na(dpID)){
    dpID <- dpID
    }else{
    if(product == "plant_diversity") dpID <- "DP1.10058.001"
    if(product == "litterfall") dpID <- "DP1.10033.001"
    if(product == "woody_veg_structure") dpID <- "DP1.10098.001"
    if(product == "belowground_biomass") dpID <- "DP1.10067.001"
    if(product == "herbaceous_clip") dpID <- "DP1.10023.001"
    if(product == "coarse_downed_wood") dpID <- "DP1.10014.001"
    if(product == "soil_microbe_biomass") dpID <- "DP1.10104.001"
  }
  neonUtilities::loadByProduct(dpID = dpID,
                               site = sites,
                               check.size = F) -> x
  return(x)
}


#' Convert raw NEON diversity object to longform plant cover data frame
#'
#' The diversity data from NEON comes as a list containing 2 data frames of data
#' that need to be combined, among other things. Here, we take those two data
#' frames and combine them into a longform data frame that can then be further
#' modified for analysis. Most of the unneccessary information from the raw data
#' has been removed. Column names that remain are plotID, subplotID, year,
#' taxonID, cover, scientificName, nativeStatusCode, family, and site.
#'
#' @import data.table
#' @importFrom data.table :=
#' @importFrom dtplyr lazy_dt
#' @param neon_div_object the raw diversity data downloaded using
#' neonPlantEcology::download_plant_div() or the function
#' neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".
#' @param trace_cover cover value for subplots where only occupancy was recorded
#' @param scale what level of spatial aggregation? This can be "1m", "10m", "100m", "plot",
#' which is the default, or "site".
#' @param timescale what level of temporal aggregation? can be "subannual", which
#' is only important for sites with multiple sampling bouts per year,
#' "annual" or "all" for the full time series.
#' @examples
#' data("D14")
#' lf <- npe_longform(D14)
#' @returns a data frame with each row a single observation of species cover at the
#' spatial and temporal scale chosen by the user.
#' @export
npe_longform <- function(neon_div_object,
                         trace_cover=0.5,
                         scale = "plot",
                         timescale = "annual"){
  .datatable.aware <- TRUE
  requireNamespace("data.table")
  requireNamespace("dplyr")
  requireNamespace("dtplyr")
  requireNamespace("tidyverse")
  requireNamespace("tidyr")
  requireNamespace("stringr")

  if(scale == "plot"){
    cover <- neon_div_object$div_1m2Data |>
      dtplyr::lazy_dt() |>
      dplyr::mutate(eventID = stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
      dplyr::mutate(endDate = as.Date(endDate)) |>
      dplyr::filter(divDataType == "plantSpecies") |>
      tidyr::replace_na(list(percentCover=trace_cover)) |>
      dplyr::filter(taxonID != "") |>
      dplyr::group_by(plotID, taxonID, eventID) |>
      dplyr::summarise(cover = sum(percentCover, na.rm=TRUE)/ifelse(as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) |>
      dplyr::ungroup() |>
      tibble::as_tibble()

    traces <- neon_div_object$div_10m2Data100m2Data |>
      dtplyr::lazy_dt() |>
      dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
      dplyr::mutate(endDate = as.Date(endDate)) |>
      dplyr::filter(targetTaxaPresent == "Y") |>
      dplyr::group_by(plotID, subplotID, taxonID, eventID) |>
      dplyr::summarise(cover = trace_cover,
                       scientificName = first(scientificName),
                       nativeStatusCode = first(nativeStatusCode),
                       family = first(family)) |>
      dplyr::ungroup() |>
      dplyr::filter(taxonID != "") |>
      dplyr::group_by(plotID, taxonID, eventID) |>
      dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) |>
      dplyr::ungroup() |>
      tibble::as_tibble()

    full_on_cover <- dplyr::bind_rows(cover, traces) |>
      dplyr::group_by(plotID, taxonID, eventID, nativeStatusCode, scientificName, family) |>
      dplyr::summarise(cover = sum(cover)) |>
      dplyr::ungroup() |>
      tidyr::replace_na(list(family = "Unknown")) |>
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4),
                    subplotID = "plot")

    if(timescale == "all") {
      full_on_cover <- full_on_cover |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
      year_range <- unique(full_on_cover$eventID) |>
        as.numeric() |>
        range() |>
        paste(collapse = "-")
      n_years <- length(unique(full_on_cover$eventID))
      full_on_cover <- full_on_cover |>
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                        family, site, subplotID) |>
        dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) |>
        dplyr::ungroup() |>
        dplyr::mutate(eventID = year_range)
    }
    if(timescale == "annual") {
      full_on_cover <- full_on_cover  |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),
                        sep = "\\.",remove = F)
      full_on_cover <- full_on_cover |>
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                        family, site, subplotID,eventID) |>
        dplyr::summarise(cover = max(cover, na.rm=T)) |>
        dplyr::ungroup()
    }
    return(full_on_cover)
  }

  if(scale == "site"){
    cover <- neon_div_object$div_1m2Data |>
      dtplyr::lazy_dt() |>
      dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
      dplyr::mutate(endDate = as.Date(endDate)) |>
      dplyr::filter(divDataType == "plantSpecies") |>
      tidyr::replace_na(list(percentCover=trace_cover)) |>
      dplyr::filter(taxonID != "") |>
      dplyr::group_by(plotID, taxonID, eventID) |>
      dplyr::summarise(cover = sum(percentCover, na.rm=TRUE)/ifelse(as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) |>
      dplyr::ungroup() |>
      tibble::as_tibble()

    traces <- neon_div_object$div_10m2Data100m2Data |>
      dtplyr::lazy_dt() |>
      dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
      dplyr::mutate(endDate = as.Date(endDate)) |>
      dplyr::filter(targetTaxaPresent == "Y") |>
      dplyr::group_by(plotID, subplotID, taxonID, eventID) |>
      dplyr::summarise(cover = trace_cover,
                       scientificName = first(scientificName),
                       nativeStatusCode = first(nativeStatusCode),
                       family = first(family)) |>
      dplyr::ungroup() |>
      dplyr::filter(taxonID != "") |>
      dplyr::group_by(plotID, taxonID, eventID) |>
      dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) |>
      dplyr::ungroup() |>
      tibble::as_tibble()

    n_plots <- length(unique(cover$plotID))

    full_on_cover <- dplyr::bind_rows(cover, traces) |>
      dtplyr::lazy_dt() |>
      dplyr::group_by(plotID, taxonID, eventID, nativeStatusCode, scientificName, family) |>
      dplyr::summarise(cover = sum(cover)) |>
      dplyr::ungroup() |>
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) |>
      dplyr::group_by(site, taxonID, eventID, nativeStatusCode, scientificName, family) |>
      dplyr::summarise(cover = sum(cover)/n_plots) |>
      dplyr::mutate(subplotID = "site",
                    plotID = "site") |>
      dplyr::ungroup() |>
      tidyr::replace_na(list(family = "Unknown")) |>
      tibble::as_tibble()

    if(timescale == "all") {
      full_on_cover <- full_on_cover  |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
      year_range <- unique(full_on_cover$eventID) |>
        as.numeric() |>
        range() |>
        paste(collapse = "-")
      n_years <- length(unique(full_on_cover$eventID))
      full_on_cover <- full_on_cover |>
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                        family, site, subplotID) |>
        dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) |>
        dplyr::ungroup() |>
        dplyr::mutate(eventID = year_range)
    }
    if(timescale == "annual") {
      full_on_cover <- full_on_cover  |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
      full_on_cover <- full_on_cover |>
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                        family, site, subplotID,eventID) |>
        dplyr::summarise(cover = max(cover, na.rm=T)) |>
        dplyr::ungroup()
    }
    return(full_on_cover)
  }

  # cover 8 ===========
  cover8 <- neon_div_object$div_1m2Data |>
    dtplyr::lazy_dt() |>
    dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
    dplyr::filter(divDataType == "plantSpecies") |>
    tidyr::replace_na(list(percentCover=trace_cover)) |>
    dplyr::select(plotID, subplotID, taxonID, eventID, cover = percentCover,
                  nativeStatusCode, scientificName, family) |>
    dplyr::filter(taxonID != "") |>
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) |>
    tidyr::replace_na(list(family = "Unknown")) |>
    tibble::as_tibble()


  # 10m2,100m2 are given 0.5 (we can change later)
  # unique(x$div_10m2Data100m2Data$subplotID) # there are 12 subplots

  # traces8 (10m2) ==============
  traces8 <- neon_div_object$div_10m2Data100m2Data |>
    dtplyr::lazy_dt() |>
    dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
    dplyr::group_by(plotID, subplotID, taxonID, eventID, scientificName,
                    nativeStatusCode, family) |>
    dplyr::summarise(cover = trace_cover) |>
    dplyr::ungroup() |>
    dplyr::filter(taxonID != "",
                  subplotID != "31", # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
                  subplotID != "32",
                  subplotID != "40",
                  subplotID != "41")  |>
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) |>
    tidyr::replace_na(list(family = "Unknown")) |>
    tibble::as_tibble()

  # traces100s ========
  traces100s <- neon_div_object$div_10m2Data100m2Data |>
    dtplyr::lazy_dt() |>
    dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
    dplyr::filter(targetTaxaPresent == "Y") |>
    dplyr::group_by(plotID, subplotID, taxonID, eventID, scientificName,
                    nativeStatusCode, family) |>
    dplyr::summarise(cover = trace_cover) |>
    dplyr::ungroup() |>
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) |>
    dplyr::filter(taxonID != "",
                  subplotID == "31"| # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
                    subplotID == "32"|
                    subplotID == "40"|
                    subplotID == "41") |>
    tidyr::replace_na(list(family = "Unknown")) |>
    tibble::as_tibble()

  # aggregating at different spatial scales ------------------------------------
  cover8_1m2 <- cover8 |>
    dplyr::group_by(plotID, subplotID, taxonID, eventID, nativeStatusCode, scientificName, family) |>
    dplyr::summarise(cover = sum(cover)) |>
    dplyr::ungroup() |>
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4))

  cover8_1m2_10m2 <- dplyr::bind_rows(cover8, traces8) |>
    dplyr::group_by(plotID,subplotID, taxonID, eventID, nativeStatusCode, scientificName, family) |>
    dplyr::summarise(cover = sum(cover)) |>
    dplyr::ungroup() |>
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4))

  cover4 <- cover8_1m2_10m2 |>
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1,2)) |>
    dplyr::bind_rows(traces100s) |> # adding in the 100m2 subplots
    dplyr::group_by(plotID, subplotID, eventID, taxonID) |>
    dplyr::summarise(cover = sum(cover), # this is summing together repeats from the rbinding
                     scientificName = first(scientificName),
                     nativeStatusCode = first(nativeStatusCode),
                     family = first(family),
                     site = first(site)) |>
    dplyr::ungroup()


  if(scale == "1m") full_on_cover <- cover8_1m2
  if(scale == "10m") full_on_cover <- cover8_1m2_10m2
  if(scale == "100m") full_on_cover <- cover4

  if(timescale == "all") {
    full_on_cover <- full_on_cover |>
      tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
    year_range <- unique(full_on_cover$eventID) |>
      as.numeric() |>
      range() |>
      paste(collapse = "-")
    n_years <- length(unique(full_on_cover$eventID))
    full_on_cover <- full_on_cover |>
      dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                      family, site, subplotID) |>
      dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) |>
      dplyr::ungroup() |>
      dplyr::mutate(eventID = year_range)
  }
  if(timescale == "annual") {
    full_on_cover <- full_on_cover |>
      tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
    full_on_cover <- full_on_cover |>
      dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName,
                      family, site, subplotID,eventID) |>
      dplyr::summarise(cover = max(cover, na.rm=T)) |>
      dplyr::ungroup()
  }

  return(full_on_cover)
}

#' Get ground cover and other variables
#'
#' @import data.table
#' @importFrom data.table :=
#' @importFrom dtplyr lazy_dt
#' @param neon_div_object  the raw diversity data downloaded using
#' neonPlantEcology::download_plant_div() or the function
#' neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".
#' @param scale the spatial scale of aggregation. Can be "1m", "10m", "100m",
#' "plot" or "site". default is "plot".
#' @param timescale The temporal scale of aggregation. Can be "all", "annual" or
#' "subannual" in the case of multiple sampling bouts per year. Defaults to "annual".
#' @examples
#' data("D14")
#' heights <- npe_groundcover(D14)
#' @returns a data frame with each row a single observation of ground cover at the
#' spatial and temporal scale chosen by the user.
#' @export
npe_groundcover <- function(neon_div_object,
                            scale = "plot",
                            timescale = "annual"){

  .datatable.aware <- TRUE
  requireNamespace("data.table")
  requireNamespace("dplyr")
  requireNamespace("dtplyr")
  requireNamespace("tidyverse")
  requireNamespace("tidyr")
  requireNamespace("stringr")

  if(scale == "plot"){
    full_on_cover <- neon_div_object$div_1m2Data |>
      dtplyr::lazy_dt() |>
      dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
      dplyr::mutate(endDate = as.Date(endDate)) |>
      dplyr::filter(divDataType == "otherVariables") |>
      # tidyr::replace_na(list(percentCover=0.5)) |>
      dplyr::filter(otherVariables != "") |>
      dplyr::group_by(plotID, otherVariables, eventID) |>
      dplyr::summarise(cover = sum(percentCover, na.rm=TRUE)/ifelse(
        as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6)) |>
      dplyr::ungroup() |>
      tibble::as_tibble() |>
      dplyr::group_by(plotID, otherVariables, eventID) |>
      dplyr::summarise(cover = sum(cover)) |>
      dplyr::ungroup() |>
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4),
                    subplotID = "plot")

    if(timescale == "all") {
      full_on_cover <- full_on_cover |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
      year_range <- unique(full_on_cover$eventID) |>
        as.numeric() |>
        range() |>
        paste(collapse = "-")
      n_years <- length(unique(full_on_cover$eventID))
      full_on_cover <- full_on_cover |>
        dplyr::group_by(plotID, otherVariables, site, subplotID) |>
        dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) |>
        dplyr::ungroup() |>
        dplyr::mutate(eventID = year_range)
    }
    if(timescale == "annual") {
      full_on_cover <- full_on_cover  |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),
                        sep = "\\.",remove = F)
      full_on_cover <- full_on_cover |>
        dplyr::group_by(plotID, otherVariables, site, subplotID,eventID) |>
        dplyr::summarise(cover = max(cover, na.rm=T)) |>
        dplyr::ungroup()
    }
    return(full_on_cover)
  }

  if(scale == "site"){
    cover <- neon_div_object$div_1m2Data |>
      dtplyr::lazy_dt() |>
      dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
      dplyr::mutate(endDate = as.Date(endDate)) |>
      dplyr::filter(divDataType == "otherVariables") |>
      # tidyr::replace_na(list(percentCover=trace_cover)) |>
      dplyr::filter(otherVariables != "") |>
      dplyr::group_by(plotID, otherVariables, eventID) |>
      dplyr::summarise(cover = sum(percentCover, na.rm=TRUE)/ifelse(
        as.numeric(stringr::str_sub(eventID,8,11))< 2019, 8,6)) |>
      dplyr::ungroup() |>
      tibble::as_tibble()

    n_plots <- length(unique(cover$plotID))

    full_on_cover <- cover |>
      dtplyr::lazy_dt() |>
      dplyr::group_by(plotID, otherVariables, eventID) |>
      dplyr::summarise(cover = sum(cover)) |>
      dplyr::ungroup() |>
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) |>
      dplyr::group_by(site, otherVariables, eventID) |>
      dplyr::summarise(cover = sum(cover)/n_plots) |>
      dplyr::mutate(subplotID = "site",
                    plotID = "site") |>
      dplyr::ungroup() |>
      tibble::as_tibble()

    if(timescale == "all") {
      full_on_cover <- full_on_cover  |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
      year_range <- unique(full_on_cover$eventID) |>
        as.numeric() |>
        range() |>
        paste(collapse = "-")
      n_years <- length(unique(full_on_cover$eventID))
      full_on_cover <- full_on_cover |>
        dplyr::group_by(plotID, otherVariables, site, subplotID) |>
        dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) |>
        dplyr::ungroup() |>
        dplyr::mutate(eventID = year_range)
    }
    if(timescale == "annual") {
      full_on_cover <- full_on_cover  |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
      full_on_cover <- full_on_cover |>
        dplyr::group_by(plotID, otherVariables, site, subplotID,eventID) |>
        dplyr::summarise(cover = max(cover, na.rm=T)) |>
        dplyr::ungroup()
    }
    return(full_on_cover)
  }

  # cover 8 ===========
  cover8 <- neon_div_object$div_1m2Data |>
    dtplyr::lazy_dt() |>
    dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
    dplyr::filter(divDataType == "otherVariables") |>
    # tidyr::replace_na(list(percentCover=trace_cover)) |>
    dplyr::select(plotID, subplotID, otherVariables, eventID, cover = percentCover) |>
    dplyr::filter(otherVariables != "") |>
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) |>
    tibble::as_tibble()

    # aggregating at different spatial scales ------------------------------------
  cover8_1m2 <- cover8 |>
    dplyr::group_by(plotID, subplotID, otherVariables, eventID) |>
    dplyr::summarise(cover = sum(cover)) |>
    dplyr::ungroup() |>
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4))

  cover4 <- cover8_1m2 |>
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1,2)) |>
    dplyr::group_by(site, plotID, subplotID, eventID, otherVariables) |>
    dplyr::summarise(cover = sum(cover)) |> # this is summing together repeats from the rbinding
    dplyr::ungroup()


  if(scale == "1m") full_on_cover <- cover8_1m2
  if(scale == "10m") full_on_cover <- cover8_1m2 # it's the same
  if(scale == "100m") full_on_cover <- cover4

  if(timescale == "all") {
    full_on_cover <- full_on_cover |>
      tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
    year_range <- unique(full_on_cover$eventID) |>
      as.numeric() |>
      range() |>
      paste(collapse = "-")
    n_years <- length(unique(full_on_cover$eventID))
    full_on_cover <- full_on_cover |>
      dplyr::group_by(plotID, otherVariables, site, subplotID) |>
      dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) |>
      dplyr::ungroup() |>
      dplyr::mutate(eventID = year_range)
  }
  if(timescale == "annual") {
    full_on_cover <- full_on_cover |>
      tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
    full_on_cover <- full_on_cover |>
      dplyr::group_by(plotID, otherVariables, site, subplotID,eventID) |>
      dplyr::summarise(cover = max(cover, na.rm=T)) |>
      dplyr::ungroup()
  }

  return(full_on_cover)
}

#' Get heights
#'
#' @import data.table
#' @importFrom data.table :=
#' @importFrom dtplyr lazy_dt
#' @param neon_div_object  the raw diversity data downloaded using
#' neonPlantEcology::download_plant_div() or the function
#' neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".
#' @param scale the spatial scale of aggregation. Can be "1m", "10m", "100m",
#' "plot" or "site". default is "plot".
#' @param timescale The temporal scale of aggregation. Can be "all", "annual" or
#' "subannual" in the case of multiple sampling bouts per year. Defaults to "annual".
#' @returns a data frame with each row a single observation of species height at the
#' spatial and temporal scale chosen by the user.
#' @export
npe_heights <- function(neon_div_object,
                            scale = "plot",
                            timescale = "annual"){

  .datatable.aware <- TRUE
  requireNamespace("data.table")
  requireNamespace("dplyr")
  requireNamespace("dtplyr")
  requireNamespace("tidyverse")
  requireNamespace("tidyr")
  requireNamespace("stringr")

  if(scale == "plot"){
    full_on_height <- neon_div_object$div_1m2Data |>
      dtplyr::lazy_dt() |>
      dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
      dplyr::mutate(endDate = as.Date(endDate)) |>
      dplyr::filter(divDataType == "plantSpecies") |>
      # tidyr::replace_na(list(percentCover=0.5)) |>
      dplyr::filter(taxonID != "") |>
      dplyr::group_by(plotID, taxonID, eventID) |>
      dplyr::summarise(
        height = mean(heightPlantSpecies, na.rm=TRUE),
        all_na = all(is.na(heightPlantSpecies))) |>
      dplyr::ungroup()  |>
      dplyr::group_by(plotID, taxonID, eventID, all_na) |>
      dplyr::summarise(height = sum(height)) |>
      dplyr::ungroup() |>
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4),
                    subplotID = "plot")  |>
      dplyr::filter(!all_na) |>
      dplyr::select(-all_na) |>
      tibble::as_tibble()

    if(timescale == "all") {
      full_on_height <- full_on_height |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
      year_range <- unique(full_on_height$eventID) |>
        as.numeric() |>
        range() |>
        paste(collapse = "-")
      n_years <- length(unique(full_on_height$eventID))
      full_on_height <- full_on_height |>
        dplyr::group_by(plotID, taxonID, site, subplotID) |>
        dplyr::summarise(height = mean(height, na.rm=T)) |>
        dplyr::ungroup() |>
        dplyr::mutate(eventID = year_range)
    }
    if(timescale == "annual") {
      full_on_height <- full_on_height  |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),
                        sep = "\\.",remove = F)
      full_on_height <- full_on_height |>
        dplyr::group_by(plotID, taxonID, site, subplotID,eventID) |>
        dplyr::summarise(height = max(height, na.rm=T)) |>
        dplyr::ungroup()
    }
    return(full_on_height)
  }

  if(scale == "site"){
    cover <- neon_div_object$div_1m2Data |>
      dtplyr::lazy_dt() |>
      dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
      dplyr::mutate(endDate = as.Date(endDate)) |>
      dplyr::filter(divDataType == "plantSpecies") |>
      dplyr::filter(taxonID != "") |>
      dplyr::group_by(plotID, taxonID, eventID) |>
      dplyr::summarise(height = mean(heightPlantSpecies, na.rm=TRUE)) |>
      dplyr::ungroup() |>
      tibble::as_tibble()

    full_on_height <- cover |>
      dtplyr::lazy_dt() |>
      dplyr::group_by(plotID, taxonID, eventID) |>
      dplyr::summarise(height = mean(height)) |>
      dplyr::ungroup() |>
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) |>
      dplyr::group_by(site, taxonID, eventID) |>
      dplyr::summarise(height = mean(height, na.rm=T)) |>
      dplyr::mutate(subplotID = "site",
                    plotID = "site") |>
      dplyr::ungroup() |>
      tibble::as_tibble()

    if(timescale == "all") {
      full_on_height <- full_on_height  |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
      year_range <- unique(full_on_height$eventID) |>
        as.numeric() |>
        range() |>
        paste(collapse = "-")
      full_on_height <- full_on_height |>
        dplyr::group_by(plotID, taxonID, site, subplotID) |>
        dplyr::summarise(cover = mean(cover, na.rm=T)) |>
        dplyr::ungroup() |>
        dplyr::mutate(eventID = year_range)
    }
    if(timescale == "annual") {
      full_on_height <- full_on_height  |>
        tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
      full_on_height <- full_on_height |>
        dplyr::group_by(plotID, taxonID, site, subplotID,eventID) |>
        dplyr::summarise(height = max(height, na.rm=T)) |>
        dplyr::ungroup()
    }
    return(full_on_height)
  }

  # cover 8 ===========
  cover8_1m2 <- neon_div_object$div_1m2Data |>
    dtplyr::lazy_dt() |>
    dplyr::mutate(eventID =stringr::str_remove_all(eventID, "\\_\\d{3}")) |>
    dplyr::filter(divDataType == "plantSpecies") |>
    dplyr::select(plotID, subplotID, taxonID, eventID, height = heightPlantSpecies) |>
    dplyr::filter(taxonID != "") |>
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) |>
    dplyr::group_by(plotID, subplotID, taxonID, eventID) |>
    dplyr::summarise(height = mean(height)) |>
    dplyr::ungroup() |>
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) |>
    tibble::as_tibble()

  cover4 <- cover8_1m2 |>
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1,2)) |>
    dplyr::group_by(site, plotID, subplotID, eventID, taxonID) |>
    dplyr::summarise(height = mean(height)) |> # this is summing together repeats from the rbinding
    dplyr::ungroup()


  if(scale == "1m") full_on_height <- cover8_1m2
  if(scale == "10m") full_on_height <- cover8_1m2 # it's the same
  if(scale == "100m") full_on_height <- cover4

  if(timescale == "all") {
    full_on_height <- full_on_height |>
      tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
    year_range <- unique(full_on_height$eventID) |>
      as.numeric() |>
      range() |>
      paste(collapse = "-")
    full_on_height <- full_on_height |>
      dplyr::group_by(plotID, taxonID, site, subplotID) |>
      dplyr::summarise(height = mean(height, na.rm=T)) |>
      dplyr::ungroup() |>
      dplyr::mutate(eventID = year_range)
  }
  if(timescale == "annual") {
    full_on_height <- full_on_height |>
      tidyr::separate(eventID, into = c("site_plot", "bout", "eventID"),sep = "\\.",remove = F)
    full_on_height <- full_on_height |>
      dplyr::group_by(plotID, taxonID, site, subplotID,eventID) |>
      dplyr::summarise(height = max(height, na.rm=T)) |>
      dplyr::ungroup()
  }

  return(full_on_height)
}
#' Create a species abundance or occurrence matrix
#'
#' npe_community_matrix creates a wide matrix of species cover or binary (presence/absence)
#' values with the plot/subplot/year as rownames. This is useful for the vegan
#' package, hence the name.
#' @param x Input object. See input argument help for more details.
#' @param neon_div_object the raw diversity data downloaded using
#' neonPlantEcology::download_plant_div() or the function
#' neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".
#' @param scale what level of aggregation? This can be "1m", "10m", "100m",
#' "plot", which is the default, or "site".
#' @param timescale what temporal resolution? can be "subannual", which is really
#' only applicable at sites where there are multiple bouts per year, "annual" or
#' "all", which dissolves together the entire time series.
#' @param trace_cover cover value for subplots where only occupancy was recorded
#' @param binary should the matrix be converted from percent cover to binary?
#' @param input by default, longform dataframe is calculated from the diversity
#' object and then converted to a community matrix, set this option to "lf"
#' to use a longform data frame that was created separately (and perhaps modified).
#' Another option is input = "divStack", which is using the output from the
#' divStack function in the neonPlants package. Using a premade longform data
#' frame or a divStack output will use the spatial and temporal scale of that
#' input data separately
#' @examples
#' data("D14")
#' heights <- npe_heights(D14)
#'
#' @returns a data frame with each row a site aggregated at the spatial and
#' temporal scales chosen by the user. Each column is a single species, and cell
#' values can be either cover (a value between 0 and 100) or occurrence (1 or 0)
#' @export
npe_community_matrix <- function(x,
                   scale = "plot",
                   trace_cover = 0.5,
                   timescale = "annual",
                   input = "neon_div_object",
                   binary=FALSE) {
  requireNamespace("tidyr")
  requireNamespace("dplyr")
  requireNamespace("tibble")

  if(input == "divStack"){

    if(scale == "plot"){
      cm <- x |>
        dplyr::select(plotID, subplotID, eventID, taxonID, targetTaxaPresent) |>
        dplyr::mutate(subplotID = "plot")
    }else{
      cm <- x |>
        dplyr::select(plotID, subplotID, eventID, taxonID, targetTaxaPresent)
    }

    if(timescale == "annual"){
      cmt <- cm |>
        dtplyr::lazy_dt() |>
        dplyr::mutate(bout = stringr::str_extract(eventID,".[\\d].") |>
                        stringr::str_remove_all("\\."),
               eventID = stringr::str_extract(eventID, "\\d{4}")) |>
        dplyr::group_by(plotID, subplotID, taxonID, eventID) |>
        dplyr::mutate(present = ifelse(targetTaxaPresent == "Y", 1, 0)) |>
        dplyr::summarise(present = max(present)) |>
        dplyr::ungroup() |>
        tibble::as_tibble()
    }
    if(timescale == "subannual"){
      cmt <- cm |>
        dtplyr::lazy_dt() |>
        dplyr::group_by(plotID, subplotID, taxonID, eventID) |>
        dplyr::mutate(present = ifelse(targetTaxaPresent == "Y", 1, 0)) |>
        dplyr::summarise(present = max(present)) |>
        dplyr::ungroup() |>
        tibble::as_tibble()
    }
    if(timescale == "all"){
      cmt <- cm |>
        dtplyr::lazy_dt() |>
        dplyr::group_by(plotID, subplotID, taxonID) |>
        dplyr::mutate(present = ifelse(targetTaxaPresent == "Y", 1, 0)) |>
        dplyr::summarise(present = max(present)) |>
        dplyr::ungroup() |>
        tibble::as_tibble()
    }

    return(
      cmt |>
        dplyr::transmute(row = paste(plotID, subplotID, eventID, sep = "_"),
                  taxonID = taxonID,
                  present=present) |>
        tidyr::pivot_wider(values_from = present, names_from = taxonID,
                           values_fill = 0,
                  values_fn = function(x)sum(x)) |>
        tibble::column_to_rownames("row")
    )
  }

  if(input == "neon_div_object"){
    longform_df <- x |>
     npe_longform(scale = scale, trace_cover = trace_cover,
                  timescale = timescale)
  }

  if(input == "lf"){longform_df <- x}

  if(!binary){
    return(
      longform_df |>
        dplyr::mutate(p_sp_y = paste(plotID, subplotID, eventID, sep = "_")) |>
        dplyr::select(p_sp_y, taxonID, cover) |>
        na.omit() |> # not sure how, but there are some NA's where they shouldn't be
        tidyr::pivot_wider(id_cols = p_sp_y,
                    names_from = taxonID,
                    values_from = cover,
                    values_fill = list(cover=0)) |>
        tibble::column_to_rownames("p_sp_y") |>
        as.data.frame()
    )
  }else{
      bin <- longform_df |>
        dplyr::mutate(p_sp_y = paste(plotID, subplotID, eventID, sep = "_")) |>
        dplyr::select(p_sp_y, taxonID, cover) |>
        na.omit() |> # not sure how, but there are some NA's where they shouldn't be
        tidyr::pivot_wider(id_cols = p_sp_y,
                    names_from = taxonID,
                    values_from = cover,
                    values_fill = list(cover=0)) |>
        tibble::column_to_rownames("p_sp_y") |>
        dplyr::mutate_all(function(x) ifelse(x>0,1,0)) |>
        as.data.frame()

      return(bin)
  }
}


#' Get plant biodiversity information for NEON plots
#'
#' npe_diversity_info calculates various biodiversity and cover indexes at the
#' plot or subplot scale at each timestep for each plot. Outputs a data frame
#' with number of species, percent cover, relative percent cover (relative to the cover of the other plants), and shannon
#' diversity, for natives, exotics and all species. Also calculates all of these
#' metrics for the families and/or species of your choice.
#'
#' @param neon_div_object the raw vegan::diversity data downloaded using
#' neonPlantEcology::download_plant_div() or #' the function
#' neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".
#' @param scale what level of aggregation? This can be "1m", "10m", "100m",
#' "plot" or "site". "plot" is the default.
#' @param trace_cover cover value for subplots where only occupancy was recorded
#' @param timescale by default npe_diversity_info groups everything by year.
#' The user may set this argument to "all" to have the function aggregate the years
#' together and then calculate diversity and cover indexes, or "subannual" for bout-level.
#' @param betadiversity If evaluating at the plot or site level, should beta
#' diversity (turnover and nestedness) be calculated. If scale = plot, it will
#' calculate betadiversity within each plot, using the combined species
#' presences within the 1 and 10 m subplots, and so it's calcuated from 8 subplots
#' before 2020, 6 after. if scale = site, it calculates the betadiversity between
#' plots.
#' @param families Which specific families should the metrics be calculated for?
#' This can be a concatenated vector if the user want more than one family.
#' @examples
#' # x <- download_plant_div("SRER")
#' # plot_level <- neonPlantEcology::npe_diversity_info(neon_div_object = x, scale = "plot")
#' @returns a data frame of higher-level summary information. Number of species,
#' Shannon-Weaver alpha diversity, cover, relative cover, for all species together
#' and grouped by nativeStatusCode.
#' @export
npe_diversity_info <- function(neon_div_object,
                               scale = "plot",
                               trace_cover = 0.5,
                               timescale = "annual",
                               betadiversity = FALSE,
                               families = NA) {
  requireNamespace("tidyr")
  requireNamespace("dplyr")
  requireNamespace('vegan')
  requireNamespace("stringr")
  # Data wrangling =============================================================

  full_on_cover <- npe_longform(neon_div_object,
                                      scale = scale,
                                      trace_cover = trace_cover,
                                      timescale = timescale)
  template <- full_on_cover |>
    dplyr::select(site, plotID, subplotID, eventID)

  # Betadiversity ===================
  if(betadiversity == TRUE & scale == "plot"){

     ten_m <- npe_longform(neon_div_object,
                           scale = "10m",
                           timescale = timescale) |>
      dplyr::group_by(site, plotID, subplotID,taxonID, eventID) |>
      dplyr::summarise(cover = sum(cover, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::group_by(site, plotID, subplotID, eventID) |>
      tidyr::spread(taxonID, cover, fill=0) |>
      dplyr::ungroup() |>
      dplyr::select(-subplotID)

     bd<- data.frame(turnover = NA, nestedness = NA, eventID = NA, plotID = NA, site = NA, subplotID = NA)

     counter <- 1
     for(i in unique(ten_m$eventID)){
       for(j in unique(ten_m$plotID)){

         if(nrow(ten_m |>
                 dplyr::filter(eventID == i, plotID == j))>0){
          out <- ten_m |>
           dplyr::filter(eventID == i, plotID == j) |>
           dplyr::select(-eventID, -site, -plotID) |>
           vegan::nestedbetajac()

          bd[counter, 1] <- out[1] |> unname()
          bd[counter, 2] <- out[2] |> unname()
          bd[counter, 3] <- i
          bd[counter, 4] <- j
          bd[counter, 5] <- stringr::str_sub(j, 1,4)
          bd[counter, 6] <- "plot"

          counter <- counter+1}
       }
     }


  }

  if(betadiversity == TRUE & scale == "site"){

    plot_scale <- npe_longform(neon_div_object,
                                scale = "plot",
                                timescale = timescale) |>
      dplyr::group_by(site, plotID,taxonID, eventID) |>
      dplyr::summarise(cover = sum(cover, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::group_by(site, plotID, eventID) |>
      tidyr::spread(taxonID, cover, fill=0) |>
      dplyr::ungroup() |>
      dplyr::select(-plotID)

    bd<- data.frame(turnover = NA, nestedness = NA, eventID = NA, site = NA, plotID = NA, subplotID = NA)

    counter <- 1
    for(i in unique(plot_scale$eventID)){
      for(j in unique(plot_scale$site)){

        if(nrow(plot_scale |>
                dplyr::filter(eventID == i, site == j))>0){
          out <- plot_scale |>
            dplyr::filter(eventID == i, site == j) |>
            dplyr::select(-eventID, -site) |>
            vegan::nestedbetajac()

          bd[counter, 1] <- out[1] |> unname()
          bd[counter, 2] <- out[2] |> unname()
          bd[counter, 3] <- i
          bd[counter, 4] <- j
          bd[counter, 5] <- "site"
          bd[counter, 6] <- "site"

          counter <- counter+1}
      }
    }

  }

  # Native vs Invasive cover ===================================================

  n_i <- full_on_cover |>
    dplyr::filter(nativeStatusCode %in% c("I", "N", "UNK")) |>
    dplyr::group_by(site, plotID, subplotID, eventID) |>
    dplyr::mutate(total_cover = sum(cover)) |>
    dplyr::ungroup() |>
    dplyr::group_by(site, plotID, subplotID,eventID, nativeStatusCode) |>
    dplyr::summarise(cover = sum(cover),
              total_cover = first(total_cover)) |>
    dplyr::ungroup() |>
    dplyr::mutate(rel_cover = cover/total_cover) |>
    dplyr::ungroup()

  lut_nsc <-c("cover_native", "cover_exotic", "cover_unknown")
  names(lut_nsc) <-  c("N", "I", "UNK")

  n_i_cover <- n_i |>
    dplyr::select(site, plotID, subplotID,eventID, nativeStatusCode, cover) |>
    dplyr::mutate(nativeStatusCode = lut_nsc[nativeStatusCode]) |>
    tidyr::pivot_wider(names_from = nativeStatusCode,
                values_from = cover,
                values_fill = list(cover = 0))

  n_i_rel_cover <- n_i |>
    dplyr::select(site, plotID, subplotID,eventID, nativeStatusCode, rel_cover) |>
    dplyr::mutate(nativeStatusCode = lut_nsc[nativeStatusCode] |> stringr::str_c("rel_", ...= _)) |>
    tidyr::pivot_wider(names_from = nativeStatusCode,
                values_from = rel_cover,
                values_fill = list(rel_cover = 0)) |>
    dplyr::left_join(n_i_cover, by = c("site", "plotID", "subplotID","eventID"))

  if(sum(names(n_i_rel_cover) %in% "cover_exotic")==0){
    n_i_rel_cover <- n_i_rel_cover |>
      dplyr::mutate(cover_exotic = 0,
             rel_cover_exotic = 0)
  }

  # not exotic cover ===================================================
  n_e <- full_on_cover |>
    dplyr::mutate(nativeStatusCode = ifelse(nativeStatusCode !="I", "NE", "I")) |>
    dplyr::group_by(site, plotID, subplotID, eventID) |>
    dplyr::mutate(total_cover = sum(cover)) |>
    dplyr::ungroup() |>
    dplyr::group_by(site, plotID, subplotID,eventID, nativeStatusCode) |>
    dplyr::summarise(cover = sum(cover),
                     total_cover = first(total_cover)) |>
    dplyr::ungroup() |>
    dplyr::mutate(rel_cover = cover/total_cover) |>
    dplyr::ungroup()

  lut_ne <-c("cover_notexotic", "cover_exotic")
  names(lut_ne) <-  c("NE", "I")

  n_e_cover <- n_e |>
    dplyr::select(site, plotID, subplotID,eventID, nativeStatusCode, cover) |>
    dplyr::mutate(nativeStatusCode = lut_ne[nativeStatusCode]) |>
    tidyr::pivot_wider(names_from = nativeStatusCode,
                       values_from = cover,
                       values_fill = list(cover = 0)) |>
    dplyr::select(-contains("cover_exotic"))

  n_e_rel_cover <- n_e |>
    dplyr::select(site, plotID, subplotID,eventID, nativeStatusCode, rel_cover) |>
    dplyr::mutate(nativeStatusCode = lut_ne[nativeStatusCode] |> stringr::str_c("rel_", ...= _)) |>
    tidyr::pivot_wider(names_from = nativeStatusCode,
                       values_from = rel_cover,
                       values_fill = list(rel_cover = 0)) |>
    dplyr::select(-contains("rel_cover_exotic")) |>
    dplyr::left_join(n_e_cover, by = c("site", "plotID", "subplotID","eventID"))



  # Cover by family ============================================================
  if(!is.na(families)){

    byfam <- full_on_cover |>
      dplyr::group_by(site, plotID, subplotID, eventID) |>
      dplyr::mutate(total_cover = sum(cover)) |>
      dplyr::ungroup() |>
      dplyr::group_by(site, plotID, subplotID,eventID, family) |>
      dplyr::summarise(cover = sum(cover),
                total_cover = first(total_cover)) |>
      dplyr::ungroup() |>
      dplyr::mutate(rel_cover = cover/total_cover) |>
      dplyr::ungroup() |>
      dplyr::filter(family %in% families)

    rcf<- byfam |>
      dplyr::select(site, plotID, subplotID,eventID, family, rel_cover) |>
      tidyr::pivot_wider(names_from = family,
                  names_prefix = "rel_cover_",
                  values_from = (rel_cover),
                  values_fill = list(rel_cover = 0))

    cf<- byfam |>
      dplyr::select(site, plotID, subplotID,eventID, family, cover) |>
      tidyr::pivot_wider(names_from = family,
                  names_prefix = "cover_",
                  values_from = (cover),
                  values_fill = list(cover = 0))

    nspp_byfam <- full_on_cover |>
      dplyr::filter(nativeStatusCode %in% c("I", "N", "UNK")) |>
      dplyr::group_by(site, plotID, subplotID, eventID) |>
      dplyr::mutate(total_cover = sum(cover)) |>
      dplyr::ungroup() |>
      dplyr::group_by(site, plotID, subplotID,eventID, family, nativeStatusCode) |>
      dplyr::summarise(nspp = length(unique(scientificName))) |>
      dplyr::ungroup() |>
      dplyr::filter(family %in% families) |>
      tidyr::pivot_wider(names_from = c(family, nativeStatusCode),
                  names_prefix = "nspp_",
                  values_from = (nspp),
                  values_fill = list(nspp = 0))
  }

  # by family, divided by biogeographic origin =================================
  if(!is.na(families)){
  family_stuff <- full_on_cover |>
    dplyr::filter(nativeStatusCode %in% c("I", "N", "UNK")) |>
    dplyr::group_by(site, plotID, subplotID, eventID) |>
    dplyr::mutate(total_cover = sum(cover)) |>
    dplyr::ungroup() |>
    dplyr::group_by(site, plotID, subplotID,eventID, family, nativeStatusCode) |>
    dplyr::summarise(cover = sum(cover),
              total_cover = first(total_cover)) |>
    dplyr::ungroup() |>
    dplyr::mutate(rel_cover = cover/total_cover) |>
    dplyr::ungroup() |>
    dplyr::filter(family %in% families)

  rc_ig<- family_stuff |>
    dplyr::select(site, plotID, subplotID,eventID, family, nativeStatusCode,rel_cover) |>
    dplyr::filter(nativeStatusCode == "I") |>
    tidyr::pivot_wider(names_from = family,
                names_prefix = "rc_exotic_",
                values_from = (rel_cover),
                values_fill = list(rel_cover = 0)) |>
    dplyr::select(-nativeStatusCode)

  rc_neg<- family_stuff |>
    dplyr::select(site, plotID, subplotID,eventID, family, nativeStatusCode,rel_cover) |>
    dplyr::filter(nativeStatusCode != "I") |>
    tidyr::pivot_wider(names_from = family,
                       names_prefix = "rc_notexotic_",
                       values_from = (rel_cover),
                       values_fill = list(rel_cover = 0)) |>
    dplyr::select(-nativeStatusCode)

  rc_ng<- family_stuff |>
    dplyr::select(site, plotID, subplotID,eventID, family, nativeStatusCode,rel_cover) |>
    dplyr::filter(nativeStatusCode == "N") |>
    tidyr::pivot_wider(names_from = family,
                names_prefix = "rc_native_",
                values_from = (rel_cover),
                values_fill = list(rel_cover = 0)) |>
    dplyr::select(-nativeStatusCode)

  c_ig<- family_stuff |>
    dplyr::select(site, plotID, subplotID,eventID, family, nativeStatusCode,cover) |>
    dplyr::filter(nativeStatusCode == "I") |>
    tidyr::pivot_wider(names_from = family,
                names_prefix = "cover_exotic_",
                values_from = (cover),
                values_fill = list(cover = 0)) |>
    dplyr::select(-nativeStatusCode)
  c_neg<- family_stuff |>
    dplyr::select(site, plotID, subplotID,eventID, family, nativeStatusCode,cover) |>
    dplyr::filter(nativeStatusCode != "I") |>
    tidyr::pivot_wider(names_from = family,
                       names_prefix = "cover_notexotic_",
                       values_from = (cover),
                       values_fill = list(cover = 0)) |>
    dplyr::select(-nativeStatusCode)
  c_ng<- family_stuff |>
    dplyr::select(site, plotID, subplotID,eventID, family, nativeStatusCode,cover) |>
    dplyr::filter(nativeStatusCode == "N") |>
    tidyr::pivot_wider(names_from = family,
                names_prefix = "cover_native_",
                values_from = (cover),
                values_fill = list(cover = 0)) |>
    dplyr::select(-nativeStatusCode)
  }
  # exotic  diversity and evenness ===============
  if(nrow(dplyr::filter(full_on_cover,nativeStatusCode=="I"))>0){
  vegan_friendly_div_ex <- full_on_cover |>
    dplyr::filter(nativeStatusCode %in% c("I")) |>
    dplyr::group_by(site, plotID, subplotID,taxonID, eventID) |>
    dplyr::summarise(cover = sum(cover, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::group_by(site, plotID, subplotID, eventID) |>
    tidyr::spread(taxonID, cover, fill=0) |>
    dplyr::ungroup()

  nspp_ex <- vegan_friendly_div_ex |>
    dplyr::select(site, plotID, subplotID,eventID) |>
    dplyr::mutate(shannon_exotic = vegan::diversity(vegan_friendly_div_ex |>
                                        dplyr::select(-site,
                                                      -plotID,
                                                      -subplotID,
                                                      -eventID)),
           evenness_exotic = shannon_exotic/vegan::specnumber(vegan_friendly_div_ex |>
                                                                  dplyr::select(-site, -plotID, -subplotID, -eventID)),
           nspp_exotic = vegan::specnumber(vegan_friendly_div_ex |>
                                      dplyr::select(-site,
                                                    -plotID,
                                                    -subplotID,
                                                    -eventID)))}else{
                                                      nspp_ex<-
                                                        template |>
                                                        dplyr::mutate(shannon_exotic = 0,
                                                               evenness_exotic = 0,
                                                               nspp_exotic = 0)
                                                    }

  # native diversity and evenness ===========
  vegan_friendly_div_n<- full_on_cover |>
    dplyr::filter(nativeStatusCode %in% c("N")) |>
    dplyr::group_by(site, plotID, subplotID,taxonID, eventID) |>
    dplyr::summarise(cover = sum(cover, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(taxonID = as.character(taxonID),
                  plotID = as.character(plotID)) |>
    dplyr::group_by(site, plotID, subplotID, eventID) |>
    tidyr::spread(taxonID, cover, fill=0) |>
    dplyr::ungroup()

  nspp_n <- vegan_friendly_div_n |>
    dplyr::select(site, plotID, subplotID,eventID) |>
    dplyr::mutate(shannon_native = vegan::diversity(vegan_friendly_div_n |>
                                               dplyr::select(-site,
                                                             -plotID,
                                                             -subplotID,
                                                             -eventID)),
           evenness_native = shannon_native/vegan::specnumber(vegan_friendly_div_n |>
                                                                  dplyr::select(-site, -plotID, -subplotID, -eventID)),
           nspp_native = vegan::specnumber(vegan_friendly_div_n |>
                                             dplyr::select(-site,
                                                           -plotID,
                                                           -subplotID,
                                                           -eventID)))


  # unknown diversity and evenness  ========================
  vegan_friendly_div_un<- full_on_cover |>
    dplyr::filter(nativeStatusCode %in% c("UNK")) |>
    dplyr::group_by(site, plotID, subplotID,taxonID, eventID) |>
    dplyr::summarise(cover = sum(cover, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(taxonID = as.character(taxonID),
                  plotID = as.character(plotID)) |>
    dplyr::group_by(site, plotID, subplotID, eventID) |>
    tidyr::spread(taxonID, cover, fill=0) |>
    dplyr::ungroup()

  nspp_un <- vegan_friendly_div_un |>
    dplyr::select(site, plotID, subplotID,eventID) |>
    dplyr::mutate(shannon_unknown = vegan::diversity(vegan_friendly_div_un |>
                                               dplyr::select(-site,
                                                             -plotID,
                                                             -subplotID,
                                                             -eventID)),
           evenness_unknown = shannon_unknown/vegan::specnumber(vegan_friendly_div_un |>
                                                              dplyr::select(-site, -plotID, -subplotID, -eventID)),
           nspp_unknown = vegan::specnumber(vegan_friendly_div_un |>
                                             dplyr::select(-site,
                                                           -plotID,
                                                           -subplotID,
                                                           -eventID)))
  # not exotic diversity and evenness  ====================
  vegan_friendly_div_nex <- full_on_cover |>
    dplyr::filter(nativeStatusCode != c("I")) |>
    dplyr::group_by(site, plotID, subplotID,taxonID, eventID) |>
    dplyr::summarise(cover = sum(cover, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(taxonID = as.character(taxonID),
                  plotID = as.character(plotID)) |>
    dplyr::group_by(site, plotID, subplotID, eventID) |>
    tidyr::spread(taxonID, cover, fill=0) |>
    dplyr::ungroup()

  nspp_nex <- vegan_friendly_div_nex |>
    dplyr::select(site, plotID, subplotID,eventID) |>
    dplyr::mutate(shannon_notexotic = vegan::diversity(vegan_friendly_div_nex |>
                                                dplyr::select(-site,
                                                              -plotID,
                                                              -subplotID,
                                                              -eventID)),
           evenness_notexotic = shannon_notexotic/vegan::specnumber(vegan_friendly_div_nex |>
                                                                  dplyr::select(-site, -plotID, -subplotID, -eventID)),
           nspp_notexotic = vegan::specnumber(vegan_friendly_div_nex |>
                                              dplyr::select(-site,
                                                            -plotID,
                                                            -subplotID,
                                                            -eventID)))

  # total vegan::diversity - not splitting between native status =========
  vegan_friendly_div_total <- full_on_cover |>
    dplyr::group_by(site, plotID, subplotID, taxonID, eventID) |>
    dplyr::summarise(cover = sum(cover)) |>
    dplyr::ungroup() |>
    dplyr::mutate(taxonID = as.character(taxonID),
           plotID = as.character(plotID)) |>
    dplyr::filter(nchar(as.character(taxonID))>0) |>
    dplyr::group_by(site, plotID, subplotID,eventID) |>
    tidyr::spread(taxonID, cover, fill=0) |>
    dplyr::ungroup()

  div_total <- dplyr::select(vegan_friendly_div_total, site, plotID, subplotID,eventID) |>
    dplyr::mutate(shannon_total = vegan::diversity(vegan_friendly_div_total |>
                                               dplyr::select(-site, -plotID, -subplotID, -eventID)),
           evenness_total = shannon_total/vegan::specnumber(vegan_friendly_div_total |>
                                                                dplyr::select(-site, -plotID, -subplotID, -eventID)),
           nspp_total = vegan::specnumber(vegan_friendly_div_total |>
                                             dplyr::select(-site, -plotID, -subplotID, -eventID)))

  # family diversity ===========================================================
  vegan_friendly_div_total_f <- full_on_cover |>
    dplyr::filter(!is.na(family)) |>
    dplyr::group_by(site, plotID, subplotID, family, eventID) |>
    dplyr::summarise(cover = sum(cover)) |>
    dplyr::ungroup() |>
    dplyr::group_by(site, plotID, subplotID,eventID) |>
    tidyr::spread(family, cover, fill=0) |>
    dplyr::ungroup()

  div_total_f <- dplyr::select(vegan_friendly_div_total_f, site, plotID, subplotID,eventID) |>
    dplyr::mutate(shannon_family = vegan::diversity(vegan_friendly_div_total_f |>
                                              dplyr::select(-site, -plotID, -subplotID, -eventID)),
           evenness_family = shannon_family/vegan::specnumber(vegan_friendly_div_total_f |>
                                               dplyr::select(-site, -plotID, -subplotID, -eventID)),
           nfamilies = vegan::specnumber(vegan_friendly_div_total_f |>
                                            dplyr::select(-site, -plotID, -subplotID, -eventID)))

  # joining and writing out ------------------------------------------------------
  final_table <- template |>
    dplyr::left_join(nspp_ex, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(nspp_nex, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(nspp_n, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(nspp_un, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(n_i_rel_cover, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(n_e_rel_cover, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(div_total, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(div_total_f, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::mutate(scale = scale,
                  invaded = ifelse(cover_exotic > 0, "invaded", "not_invaded"))
  if(exists("bd")){
    final_table <- final_table |>
      dplyr::left_join(bd, by = c("site", "plotID", "subplotID", "eventID"))
  }
  if(!is.na(families)){
    final_table <- final_table |>
    dplyr::left_join(rcf, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(rc_ig, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(c_ig, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(cf, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(rc_ng, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(c_ng, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(rc_neg, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(c_neg, by = c("site", "plotID", "subplotID", "eventID")) |>
    dplyr::left_join(nspp_byfam, by = c("site", "plotID", "subplotID", "eventID"))}

  # seems crazy, i know... but those NAs should all definitely be zero
  final_table <- final_table |>
    dplyr::mutate_all(list(~ replace(., is.na(.), 0))) |>
    unique() # temporary fix, for some reason it's returning repeats of each row - there's mutate somewhere where there needs to be a summarise maybe

  return(final_table)
}

#' Change the native status code for a particular taxon at a particular site
#'
#' Sometimes even though a particular species identity is not known, the end
#' user can still determine its native status. For example, maybe the taxon
#' was identified to the genus level, and the local flora confirms that all
#' plants in that genus are native at that particular site. This function
#' allows for post-hoc modification of the native status code for cases like this.
#'
#' @param df is the data frame returned by npe_longform
#' @param taxon is the taxonID column in the data frame
#' @param site is the identity of the NEON site (e.g. "JORN")
#' @param new_code is the NativeStatusCode value to change to
#' @returns a data frame
#' @examples
#'
#' data("D14")
#' # convert to longform cover
#' lf_div <- npe_longform(D14)
#'
#' # change all of the unknown Abutilon spp to native
#' modified_lf_div <- npe_change_native_status(lf_div, "ABUTI", "JORN", "N")
#'
#' @export
npe_change_native_status <- function(df, taxon, site, new_code){
  requireNamespace('dplyr')
  return(
    df |>
      dplyr::mutate(nativeStatusCode = replace(nativeStatusCode,
                                        taxonID == taxon & site == site,
                                        new_code))
  )
}



#'Download and join spatial information to a neonPlantEcology output data frame
#'
#'@param df a neonPlantEcology-produced data frame
#'@param type what type of ancillary data structure you want joined. Can be
#'"spatial", which will turn the data frame into an sf data frame, or "latlong",
#'which will add the latitudes and longitudes and other ancillary data as
#'columns only.
#'@param spatial_only set to TRUE if you only want the coordinates and none of
#'the ancillary variables.
#'@param dest_dir where to download the files
#'@param input to what kind of neonPlantEcology product are you appending? Can
#'be "community_matrix", "longform_cover", or "summary_info".
#'@returns a data frame
#'
#'@export
npe_plot_centroids <- function(df,
                               dest_dir = file.path(getwd(), "tmp"),
                               type = "latlong",
                               spatial_only = TRUE,
                               input = "community_matrix"){
  requireNamespace("stringr")
  requireNamespace("sf")
  requireNamespace("tibble")
  requireNamespace("dplyr")
  shp_file <- file.path(dest_dir, "All_NEON_TOS_Plots_V8/All_NEON_TOS_Plot_Centroids_V8.shp")
  if(!file.exists(shp_file)){
    url <- "https://www.neonscience.org/sites/default/files/All_NEON_TOS_Plots_V8.zip"
    file <- stringr::str_split(url, "/", simplify = T)[length(stringr::str_split(url, "/", simplify = T))]
    dir.create(dest_dir, recursive = T)
    download.file(url=url, destfile = file.path(dest_dir, file))
    unzip(zipfile = file.path(dest_dir,file),
          exdir = dest_dir)
  }
  if(type == "spatial") neon_plots <- sf::st_read(shp_file)
  if(type == "latlong") neon_plots <- sf::st_read(stringr::str_replace(shp_file, ".shp", ".csv"))

  if(input == "community_matrix") outdf <- df |>
      tibble::rownames_to_column("plot_info") |>
      dplyr::mutate(plotID = stringr::str_sub(plot_info,1,8)) |>
      dplyr::left_join(neon_plots, by = "plotID")
  if(input == "longform_cover") outdf <- df |>
      dplyr::left_join(neon_plots, by = "plotID")
  if(input == "summary_info") outdf <- df |>
      dplyr::left_join(neon_plots, by = "plotID")

  if(spatial_only && type == "spatial") outdf <- dplyr::select(outdf, plotID)
  if(spatial_only && type == "latlong") outdf <- dplyr::select(outdf, plotID, latitude, longitude)
  return(outdf)
}


#' Get plot information from a community matrix
#'
#' The npe_community_matrix() function is designed to work with the vegan
#' package, and one of the requirements of vegan functions is that there are
#' only numeric columns in community matrices. Therefore, all of the metatdata
#' is collapsed into the rownames. This function allows you to extract that
#' very basic metadata back out to a more easily interpretable data frame.
#'
#'@param comm the community matrix object created by npe_community_matrix()
#'@returns a data frame
#'@examples
#'data("D14")
#'npe_community_matrix(D14) |> npe_plot_info()
#'@export
npe_plot_info <- function(comm){
  requireNamespace("dplyr")
  requireNamespace("tibble")
  requireNamespace("tidyr")
  requireNamespace("stringr")
  return(comm |>
           tibble::rownames_to_column("rowname") |>
           tidyr::separate(rowname, into = c("site", "plot", "scale", "eventID"), remove = F) |>
           dplyr::mutate(plotID = stringr::str_c(site, "_", plot)) |>
           dplyr::select(site, plot, scale, eventID, plotID, rowname))
}

#' get site ids
#'
#' This uses the site boundary shapefile (obtainable by data('sites')) to get a
#' list of siteID codes to feed into npe_download.
#'
#' @param domain can be one or more domain codes, as a character vector, or as a number.
#' e.g. domain = c("D01", "D14"), or domain = c(3, 14), can also be a mix: domain = c(3, "D04).
#' @param type can be "Core Terrestrial" or "Relocatable Terrestrial"
#' @examples
#'
#' # if no domains or site types are specified, it returns all site codes
#' all_sites <- npe_site_ids()
#' npe_site_ids(domain = c("Northeast", "Mid-Atlantic"))
#' npe_site_ids(domain = c("D02", 15))
#'
#' @returns a vector of four letter site identification codes.
#' @export
npe_site_ids <- function(domain = NA, type = NA){
  requireNamespace("dplyr")
  requireNamespace("stringr")
  data("sites")
  # all
  if(is.na(domain[1]) & is.na(type[1])) return(unique(sites$siteID))

  # fixing common typos in domain argument
  if(any(is.numeric(domain))) domain <-
    stringr::str_c("D", stringr::str_pad(domain, width = 2, side = "left", pad = "0"))
  if(any(is.character(domain)) & any(nchar(domain) < 3)) domain <-
    stringr::str_remove_all(domain, "D") |>
    stringr::str_pad(width = 2, side = "left", pad = "0") |>
    stringr::str_c("D", pattern = _)
  if(any(is.character(domain)) & any(nchar(domain) > 3)){
    sites <- dplyr::filter(sites, domainName %in% domain)
    return(sites |> dplyr::pull(siteID) |> unique())
  }

  if(!is.na(domain[1])) sites <- dplyr::filter(sites, domainNumb %in% domain)
  if(!is.na(type[1])) sites <- dplyr::filter(sites, type %in% siteType)
  return(sites |> dplyr::pull(siteID) |> unique())
}



