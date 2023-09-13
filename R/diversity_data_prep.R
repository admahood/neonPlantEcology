
#' Plant diversity data downloader
#'
#' A function to download plant diversity data from the NEON API.
#'
#' Data comes in two separate components:
#' 1m2 subplots with cover estimates (first item of the list)
#' 10 and 100m2 subplots with presence only (second item of the list)
#' @param sites a vector of NEON site abbreviations. Defaults to "SRER"
#' @keywords download neon diversity
#'
#' @examples
#' # diversity_object <- download_plant_div(sites = "SRER")
#' @export
npe_download_plant_div <- function(sites = "SRER"){
  requireNamespace("neonUtilities")
  neonUtilities::loadByProduct(dpID = "DP1.10058.001",
                site = sites,
                check.size = F) -> x
  return(x)
}


#' Data downloader
#'
#' A wrapper function to download data from the NEON API. Some commonly used
#' products are provided as plain language options, otherwise the user
#' can enter the product ID number (dpID).
#'
#' @param sites a vector of NEON site abbreviations. Defaults to "SRER"
#' @param product a plain language vector of the data product to be downloaded.
#' Can be "plant_diversity", "litterfall", "woody_veg_structure",
#' "belowground_biomass", "herbaceous_clip", "coarse_downed_wood",
#' or "soil_microbe_biomass"
#' @param dpID if you need a data product not given as one of the product
#' options, set the data product ID here (e.g. "DP1.10023.001").
#' @keywords download neon diversity
#'
#' @examples
#' # diversity_object <- npe_download(sites = "SRER")
#' @export
npe_download <- function(sites = "SRER",
                         dpID = NA,
                         product = "plant_diversity"){
  requireNamespace("neonUtilities")

  if(!is.na(dpID)){dpID <- dpID}else{
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

############################################################
# name_cleaner tries to fix as many typos, etc as possible #
############################################################

#' Clean scientific names
#'
#' This function does some straight forward grep-ing to aggregate to clean up
#' some of the latin names
#'
#' @param lf_cover the longform cover table created by using
#' neonPlantEcology::npe_longform().
npe_name_cleaner <- function(lf_cover){
  requireNamespace('magrittr')
  requireNamespace('dplyr')
  all_sp=sort(unique(lf_cover$scientificName))

  species_names=tolower(all_sp)


  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+) \\(l.\\).+$",replacement = "\\1")
  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+-\\w+) \\(l.\\).+$",replacement = "\\1")

  species_names=gsub(x = species_names,pattern = "^(\\w+ \\w+) \\(.+\\).*$",replacement = "\\1")

  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+-\\w+) l.$",replacement = "\\1")
  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+) l\\..*$",replacement = "\\1")

  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+) \\w+ ssp.*",replacement = "\\1")
  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+) \\w+ var.*",replacement = "\\1")

  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+) \\w+\\. ssp.*",replacement = "\\1")
  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+) \\w+\\. var.*",replacement = "\\1")

  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+) \\w+$",replacement = "\\1")
  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+) \\w+\\.$",replacement = "\\1")

  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+) .*",replacement = "\\1")
  species_names=gsub(x = species_names,pattern = "(^\\w+ \\w+-\\w+) .*",replacement = "\\1")

  species_names=gsub(x = species_names,pattern = "(^\\w+ ×\\w+-\\w+).*$",replacement = "\\1")
  species_names=gsub(x = species_names,pattern = "(^\\w+ ×\\w+).*$",replacement = "\\1")


  species_names=gsub(x = species_names,pattern = "(^\\w+ )\\w+/\\w+$",replacement = "\\1sp.")
  species_names<-gsub(x = species_names,pattern = "^\\w+ \\w+/.+$",replacement = "unknown plant")

  species_names <- str_to_sentence(species_names)

  lut_sn <- all_sp
  names(lut_sn) <- species_names

  return(
    lf_cover %>%
      dplyr::mutate(scientificName = lut_sn[scientificName])
    )


}

#######################################################################
# npe_longform creates a long dataframe for cover from the neon #
# diversity object (used in all of the following functions)           #
#######################################################################

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
#' @param scale what level of aggregation? This can be "1m", "10m", "100m", "plot",
#' which is the default, or "site".
#' @param fix_unks Should the unknown codes be altered with the "unk_fixer()"
#' function? Defaults to false. This requires manual investigation and editing
#' of the unk_fixer function.
#' @examples
#' # raw_div <- npe_download_plant_div(sites = "SRER")
#' # lf_div <- npe_longform(raw_div)
#' @export
npe_longform <- function(neon_div_object,
                               trace_cover=0.5,
                               scale = "plot",
                               dissolve_years = FALSE,
                               fix_unks = FALSE){
  .datatable.aware <- TRUE
  requireNamespace("data.table")
  requireNamespace("dplyr")
  requireNamespace("dtplyr")
  requireNamespace("tidyverse")
  requireNamespace("tidyr")
  requireNamespace("stringr")
  requireNamespace("magrittr")

  if(scale == "plot"){
    cover <- neon_div_object$div_1m2Data %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(divDataType == "plantSpecies") %>%
      dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4))) %>%
      tidyr::replace_na(list(percentCover=trace_cover)) %>%
      dplyr::group_by(plotID, subplotID, taxonID, year) %>%
      # dealing with the multiple bout issue by first getting the max cover
      # per sampling effort
      dplyr::summarise(cover = max(percentCover),
                nativeStatusCode = first(nativeStatusCode),
                scientificName = first(scientificName),
                family = first(family)) %>%
      dplyr::ungroup()  %>%
      dplyr::filter(taxonID != "") %>%
      dplyr::group_by(plotID, taxonID, year) %>%
      dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(year)<2019, 8,6),
                nativeStatusCode = first(nativeStatusCode),
                scientificName = first(scientificName),
                family = first(family)) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble()

    traces <- neon_div_object$div_10m2Data100m2Data %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(targetTaxaPresent == "Y") %>%
      dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
      dplyr::group_by(plotID, subplotID, taxonID, year) %>%
      dplyr::summarise(cover = trace_cover,
                scientificName = first(scientificName),
                nativeStatusCode = first(nativeStatusCode),
                family = first(family)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(taxonID != "") %>%
      dplyr::group_by(plotID, taxonID, year) %>%
      dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(year)<2019, 12,10),
                nativeStatusCode = first(nativeStatusCode),
                scientificName = first(scientificName),
                family = first(family)) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble()

    full_on_cover <- dplyr::bind_rows(cover, traces) %>%
      dplyr::group_by(plotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
      dplyr::summarise(cover = sum(cover)) %>%
      dplyr::ungroup()%>%
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4),
             subplotID = "plot")
    if(fix_unks) full_on_cover <- full_on_cover %>%  unk_fixer()

    if(dissolve_years) {
      year_range <- unique(full_on_cover$year)%>%
        as.numeric %>%
        range %>%
        paste(collapse = "-")
      n_years <- length(unique(full_on_cover$year))
      full_on_cover <- full_on_cover %>%
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName, family, site, subplotID) %>%
        dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(year = year_range)
    }

    return(full_on_cover)
  }
  if(scale == "site"){
    cover <- neon_div_object$div_1m2Data %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(divDataType == "plantSpecies") %>%
      dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
      tidyr::replace_na(list(percentCover=trace_cover)) %>%
      dplyr::group_by(plotID, subplotID, taxonID, year) %>%
      # dealing with the multiple bout issue by first getting the max cover
      # per sampling effort
      dplyr::summarise(cover = max(percentCover),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) %>%
      dplyr::ungroup()  %>%
      dplyr::filter(taxonID != "") %>%
      dplyr::group_by(plotID, taxonID, year) %>%
      dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(year)<2019, 8,6),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble()

    traces <- neon_div_object$div_10m2Data100m2Data %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(targetTaxaPresent == "Y") %>%
      dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
      dplyr::group_by(plotID, subplotID, taxonID, year) %>%
      dplyr::summarise(cover = trace_cover,
                       scientificName = first(scientificName),
                       nativeStatusCode = first(nativeStatusCode),
                       family = first(family)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(taxonID != "") %>%
      dplyr::group_by(plotID, taxonID, year) %>%
      dplyr::summarise(cover = sum(cover, na.rm=TRUE)/ifelse(as.numeric(year)<2019, 12,10),
                       nativeStatusCode = first(nativeStatusCode),
                       scientificName = first(scientificName),
                       family = first(family)) %>%
      dplyr::ungroup() %>%
      tibble::as_tibble()

    n_plots <- length(unique(cover$plotID))

    full_on_cover <- dplyr::bind_rows(cover, traces) %>%
      dplyr::group_by(plotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
      dplyr::summarise(cover = sum(cover)) %>%
      dplyr::ungroup()%>%
      dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) %>%
      dplyr::group_by(site, taxonID, year, nativeStatusCode, scientificName, family) %>%
      dplyr::summarise(cover = sum(cover)/n_plots) %>%
      dplyr::mutate(subplotID = "site",
             plotID = "site") %>%
      dplyr::ungroup()
    if(fix_unks) full_on_cover <- full_on_cover %>%  unk_fixer()
    if(dissolve_years) {
      year_range <- unique(full_on_cover$year)%>%
        as.numeric %>%
        range %>%
        paste(collapse = "-")
      n_years <- length(unique(full_on_cover$year))
      full_on_cover <- full_on_cover %>%
        dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName, family, site, subplotID) %>%
        dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(year = year_range)
    }
    return(full_on_cover)
  }

  cover8 <- neon_div_object$div_1m2Data %>%
    dtplyr::lazy_dt() %>%
    dplyr::mutate(endDate = as.Date(endDate)) %>%
    dplyr::filter(divDataType == "plantSpecies") %>%
    dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
    # entries in the df with no values but species was there
    # i.e. someone put the sp. code and forgot to fill in the number
    # putting as trace cover value
    tidyr::replace_na(list(percentCover=trace_cover)) %>%
    dplyr::mutate(endDate = as.Date(endDate)) %>%
    dplyr::filter(divDataType == "plantSpecies") %>%
    dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
    dplyr::group_by(plotID, subplotID, taxonID, year) %>%
    # dealing with the multiple bout issue by first getting the mean cover
    # per sampling effort, without aggregating, then later we'll aggregate.
    # that way, a fall-bloomer that isn't visible in spring, for example,
    # will be given its full cover value for fall, but then a species
    # that is there for both seasons will be averaged, if that makes sense
    dplyr::summarise(cover = max(percentCover),
              nativeStatusCode = first(nativeStatusCode),
              scientificName = first(scientificName),
              family = first(family)) %>%
    dplyr::ungroup()  %>%
    dplyr::filter(taxonID != "") %>%
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) %>%
    tibble::as_tibble()


  # 10m2,100m2 are given 0.5 (we can change later)
  # unique(x$div_10m2Data100m2Data$subplotID) # there are 12 subplots

  traces8 <- neon_div_object$div_10m2Data100m2Data %>%
    dtplyr::lazy_dt() %>%
    dplyr::mutate(endDate = as.Date(endDate)) %>%
    dplyr::filter(targetTaxaPresent == "Y") %>%
    dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
    dplyr::group_by(plotID, subplotID, taxonID, year) %>%
    dplyr::summarise(cover = trace_cover,
              scientificName = first(scientificName),
              nativeStatusCode = first(nativeStatusCode),
              family = first(family)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(taxonID != "",
           subplotID != "31", # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
           subplotID != "32",
           subplotID != "40",
           subplotID != "41")  %>%
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1, 4)) %>%
    tibble::as_tibble()

  traces100s <- neon_div_object$div_10m2Data100m2Data %>%
    dtplyr::lazy_dt() %>%
    dplyr::mutate(endDate = as.Date(endDate)) %>%
    dplyr::filter(targetTaxaPresent == "Y") %>%
    dplyr::mutate(year = stringr::str_c(stringr::str_sub(endDate,1,4)))%>%
    dplyr::group_by(plotID, subplotID, taxonID, year) %>%
    dplyr::summarise(cover = trace_cover,
              scientificName = first(scientificName),
              nativeStatusCode = first(nativeStatusCode),
              family = first(family)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4)) %>%
    dplyr::filter(taxonID != "",
           subplotID == "31"| # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
             subplotID == "32"|
             subplotID == "40"|
             subplotID == "41") %>%
    tibble::as_tibble()

  # aggregating at different scales ----------------------------------------------
  cover8_1m2 <- cover8 %>%
    dplyr::group_by(plotID, subplotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
    dplyr::summarise(cover = sum(cover)) %>%
    dplyr::ungroup()%>%
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4))
  if(fix_unks) cover8_1m2 <- unk_fixer(cover8_1m2)

  cover8_1m2_10m2 <- dplyr::bind_rows(cover8, traces8) %>%
    dplyr::group_by(plotID,subplotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
    dplyr::summarise(cover = sum(cover)) %>%
    dplyr::ungroup()%>%
    dplyr::mutate(site = stringr::str_sub(plotID, 1,4))
  if(fix_unks) cover8_1m2_10m2<-cover8_1m2_10m2 %>%  unk_fixer()

  cover4 <- cover8_1m2_10m2 %>%
    dplyr::mutate(subplotID = stringr::str_sub(subplotID, 1,2)) %>%
    dplyr::bind_rows(traces100s) %>% # adding in the 100m2 subplots
    dplyr::group_by(plotID, subplotID, year, taxonID) %>%
    dplyr::summarise(cover = sum(cover), # this is summing together repeats from the rbinding
              scientificName = first(scientificName),
              nativeStatusCode = first(nativeStatusCode),
              family = first(family),
              site = first(site)) %>%
    dplyr::ungroup()
  if(fix_unks) cover4 <- cover4 %>%  unk_fixer()


  if(scale == "1m") full_on_cover <- cover8_1m2
  if(scale == "10m") full_on_cover <- cover8_1m2_10m2
  if(scale == "100m") full_on_cover <- cover4

  if(dissolve_years) {
    year_range <- unique(full_on_cover$year)%>%
      as.numeric %>%
      range %>%
      paste(collapse = "-")
    n_years <- length(unique(full_on_cover$year))
    full_on_cover <- full_on_cover %>%
      dplyr::group_by(plotID, taxonID, nativeStatusCode, scientificName, family, site, subplotID) %>%
      dplyr::summarise(cover = sum(cover, na.rm=T)/n_years) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(year = year_range)
  }

  return(full_on_cover)
}

#' Create a species abundance or occurrence matrix
#'
#' npe_community_matrix creates a wide matrix of species cover or binary (presence/absence)
#' values with the plot/subplot/year as rownames. This is useful for the vegan
#' package, hence the name.
#'
#' @param neon_div_object the raw diversity data downloaded using
#' neonPlantEcology::download_plant_div() or the function
#' neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".
#' @param scale what level of aggregation? This can be "1m", "10m", "100m", or "plot",
#' which is the default.
#' @param trace_cover cover value for subplots where only occupancy was recorded
#' @param fix_unks Should the unknown codes be altered with the "unk_fixer()"
#' function? Defaults to false. This requires manual investigation and editing
#' of the unk_fixer function.
#' @param binary should the matrix be converted from percent cover to binary?
#' @export
npe_community_matrix <- function(neon_div_object,
                   scale="plot",
                   trace_cover = 0.5,
                   fix_unks = FALSE,
                   binary=FALSE) {
  requireNamespace("tidyr")
  requireNamespace("dplyr")
  requireNamespace("tibble")
  requireNamespace("magrittr")

  if(!binary){
    return(
      neon_div_object %>%
        npe_longform(scale = scale, trace_cover = trace_cover, fix_unks = FALSE) %>%
        dplyr::mutate(p_sp_y = paste(plotID, subplotID, year, sep = "_")) %>%
        dplyr::select(p_sp_y, taxonID, cover) %>%
        na.omit() %>% # not sure how, but there are some NA's where they shouldn't be
        tidyr::pivot_wider(id_cols = p_sp_y,
                    names_from = taxonID,
                    values_from = cover,
                    values_fill = list(cover=0)) %>%
        tibble::column_to_rownames("p_sp_y") %>%
        as.data.frame()
    )
  }else{
      bin<-neon_div_object %>%
        npe_longform(scale = scale, trace_cover = trace_cover, fix_unks = FALSE) %>%
        dplyr::mutate(p_sp_y = paste(plotID, subplotID, year, sep = "_")) %>%
        dplyr::select(p_sp_y, taxonID, cover) %>%
        na.omit() %>% # not sure how, but there are some NA's where they shouldn't be
        tidyr::pivot_wider(id_cols = p_sp_y,
                    names_from = taxonID,
                    values_from = cover,
                    values_fill = list(cover=0)) %>%
        tibble::column_to_rownames("p_sp_y") %>%
        dplyr::mutate_all(function(x) ifelse(x>0,1,0)) %>%
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
#' @param scale what level of aggregation? This can be "1m", "10m", "100m", "plot",
#' which is the default, or "site".
#' @param trace_cover cover value for subplots where only occupancy was recorded
#' @param fix_unks Should the unknown codes be altered with the "unk_fixer()"
#' function? Defaults to false. This requires manual investigation and editing
#' of the unk_fixer function.
#' @param dissolve_years by default npe_diversity_info groups everything by year.
#' The user may set this argument to true to have the function aggregate the years
#' together and then calculate diversity and cover indexes.
#' @param betadiversity If evaluating at the plot or site level, should beta
#' diversity (turnover and nestedness) be calculated. If scale = plot, it will
#' calculate betadiversity within each plot, using the combined species
#' presences within the 1 and 10 m subplots, and so it's calcuated from 8 subplots
#' before 2020, 6 after. if scale = site, it calculates the betadiversity between
#' plots.
#' @param families Which specific families should the metrics be calculated for?
#' This can be a concatenated vector if the user want more than one family.
#' @param species Which specific species should the metrics be calculated for?
#' This can be a concatenated vector if the user want more than one species.
#' @examples
#' # x <- download_plant_div("SRER")
#' # plot_level <- neonPlantEcology::npe_diversity_info(neon_div_object = x, scale = "plot")
#' @export
npe_diversity_info <- function(neon_div_object,
                               scale = "plot",
                               trace_cover = 0.5,
                               fix_unks = FALSE,
                               dissolve_years = FALSE,
                               betadiversity = FALSE,
                               name_cleaner = FALSE,
                               families = NA,
                               spp = NA) {
  requireNamespace("tidyr")
  requireNamespace("dplyr")
  requireNamespace('vegan')
  requireNamespace("magrittr")
  # Data wrangling =============================================================

  full_on_cover <- npe_longform(neon_div_object,
                                      scale = scale,
                                      dissolve_years = dissolve_years,
                                      fix_unks = fix_unks)
  if(name_cleaner) full_on_cover <- name_cleaner(full_on_cover)

  template <- full_on_cover %>%
    dplyr::select(site, plotID, subplotID, year)

  # Betadiversity ===================
  if(betadiversity == TRUE & scale == "plot"){

     ten_m <- npe_longform(neon_div_object,
                                 scale = "10m",
                                 dissolve_years = dissolve_years,
                                 fix_unks = fix_unks) %>%
      dplyr::group_by(site, plotID, subplotID,taxonID, year) %>%
      dplyr::summarise(cover = sum(cover, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(site, plotID, subplotID, year) %>%
      tidyr::spread(taxonID, cover, fill=0) %>%
      dplyr::ungroup() %>%
      dplyr::select(-subplotID)

     bd<- data.frame(turnover = NA, nestedness = NA, year = NA, plotID = NA, site = NA, subplotID = NA)

     counter <- 1
     for(i in unique(ten_m$year)){
       for(j in unique(ten_m$plotID)){

         if(nrow(ten_m %>%
                 dplyr::filter(year == i, plotID == j))>0){
          out <- ten_m %>%
           dplyr::filter(year == i, plotID == j) %>%
           dplyr::select(-year, -site, -plotID) %>%
           vegan::nestedbetajac()

          bd[counter, 1] <- out[1] %>% unname
          bd[counter, 2] <- out[2] %>% unname
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
                                dissolve_years = dissolve_years,
                                fix_unks = fix_unks) %>%
      dplyr::group_by(site, plotID,taxonID, year) %>%
      dplyr::summarise(cover = sum(cover, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(site, plotID, year) %>%
      tidyr::spread(taxonID, cover, fill=0) %>%
      dplyr::ungroup() %>%
      dplyr::select(-plotID)

    bd<- data.frame(turnover = NA, nestedness = NA, year = NA, site = NA, plotID = NA, subplotID = NA)

    counter <- 1
    for(i in unique(plot_scale$year)){
      for(j in unique(plot_scale$site)){

        if(nrow(plot_scale %>%
                dplyr::filter(year == i, site == j))>0){
          out <- plot_scale %>%
            dplyr::filter(year == i, site == j) %>%
            dplyr::select(-year, -site) %>%
            vegan::nestedbetajac()

          bd[counter, 1] <- out[1] %>% unname
          bd[counter, 2] <- out[2] %>% unname
          bd[counter, 3] <- i
          bd[counter, 4] <- j
          bd[counter, 5] <- "site"
          bd[counter, 6] <- "site"

          counter <- counter+1}
      }
    }

  }

  # Native vs Invasive cover ===================================================

  n_i <- full_on_cover %>%
    dplyr::filter(nativeStatusCode %in% c("I", "N", "UNK")) %>%
    dplyr::group_by(site, plotID, subplotID, year) %>%
    dplyr::mutate(total_cover = sum(cover))%>%
    dplyr::ungroup() %>%
    dplyr::group_by(site, plotID, subplotID,year, nativeStatusCode) %>%
    dplyr::summarise(cover = sum(cover),
              total_cover = first(total_cover)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rel_cover = cover/total_cover) %>%
    dplyr::ungroup()

  lut_nsc <-c("cover_native", "cover_exotic", "cover_unknown")
  names(lut_nsc) <-  c("N", "I", "UNK")

  n_i_cover <- n_i %>%
    dplyr::select(site, plotID, subplotID,year, nativeStatusCode, cover) %>%
    dplyr::mutate(nativeStatusCode = lut_nsc[nativeStatusCode]) %>%
    tidyr::pivot_wider(names_from = nativeStatusCode,
                values_from = cover,
                values_fill = list(cover = 0))

  n_i_rel_cover <- n_i %>%
    dplyr::select(site, plotID, subplotID,year, nativeStatusCode, rel_cover) %>%
    dplyr::mutate(nativeStatusCode = lut_nsc[nativeStatusCode] %>% stringr::str_c("rel_",.)) %>%
    tidyr::pivot_wider(names_from = nativeStatusCode,
                values_from = rel_cover,
                values_fill = list(rel_cover = 0))%>%
    dplyr::left_join(n_i_cover, by = c("site", "plotID", "subplotID","year"))

  if(sum(names(n_i_rel_cover) %in% "cover_exotic")==0){
    n_i_rel_cover <- n_i_rel_cover %>%
      dplyr::mutate(cover_exotic = 0,
             rel_cover_exotic = 0)
  }

  # not exotic cover ===================================================
  n_e <- full_on_cover %>%
    dplyr::mutate(nativeStatusCode = ifelse(nativeStatusCode !="I", "NE", "I")) %>%
    dplyr::group_by(site, plotID, subplotID, year) %>%
    dplyr::mutate(total_cover = sum(cover))%>%
    dplyr::ungroup() %>%
    dplyr::group_by(site, plotID, subplotID,year, nativeStatusCode) %>%
    dplyr::summarise(cover = sum(cover),
                     total_cover = first(total_cover)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rel_cover = cover/total_cover) %>%
    dplyr::ungroup()

  lut_ne <-c("cover_notexotic", "cover_exotic")
  names(lut_ne) <-  c("NE", "I")

  n_e_cover <- n_e %>%
    dplyr::select(site, plotID, subplotID,year, nativeStatusCode, cover) %>%
    dplyr::mutate(nativeStatusCode = lut_ne[nativeStatusCode]) %>%
    tidyr::pivot_wider(names_from = nativeStatusCode,
                       values_from = cover,
                       values_fill = list(cover = 0)) %>%
    dplyr::select(-contains("cover_exotic"))

  n_e_rel_cover <- n_e %>%
    dplyr::select(site, plotID, subplotID,year, nativeStatusCode, rel_cover) %>%
    dplyr::mutate(nativeStatusCode = lut_ne[nativeStatusCode] %>% stringr::str_c("rel_",.)) %>%
    tidyr::pivot_wider(names_from = nativeStatusCode,
                       values_from = rel_cover,
                       values_fill = list(rel_cover = 0))%>%
    dplyr::select(-contains("rel_cover_exotic")) %>%
    dplyr::left_join(n_e_cover, by = c("site", "plotID", "subplotID","year"))



  # Cover by family ============================================================
  if(!is.na(families)){

    byfam <- full_on_cover%>%
      dplyr::group_by(site, plotID, subplotID, year) %>%
      dplyr::mutate(total_cover = sum(cover))%>%
      dplyr::ungroup() %>%
      dplyr::group_by(site, plotID, subplotID,year, family) %>%
      dplyr::summarise(cover = sum(cover),
                total_cover = first(total_cover)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(rel_cover = cover/total_cover) %>%
      dplyr::ungroup() %>%
      dplyr::filter(family %in% families)

    rcf<- byfam%>%
      dplyr::select(site, plotID, subplotID,year, family, rel_cover) %>%
      tidyr::pivot_wider(names_from = family,
                  names_prefix = "rel_cover_",
                  values_from = (rel_cover),
                  values_fill = list(rel_cover = 0))

    cf<- byfam%>%
      dplyr::select(site, plotID, subplotID,year, family, cover) %>%
      tidyr::pivot_wider(names_from = family,
                  names_prefix = "cover_",
                  values_from = (cover),
                  values_fill = list(cover = 0))

    nspp_byfam <- full_on_cover%>%
      dplyr::filter(nativeStatusCode %in% c("I", "N", "UNK")) %>%
      dplyr::group_by(site, plotID, subplotID, year) %>%
      dplyr::mutate(total_cover = sum(cover))%>%
      dplyr::ungroup() %>%
      dplyr::group_by(site, plotID, subplotID,year, family, nativeStatusCode) %>%
      dplyr::summarise(nspp = length(unique(scientificName))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(family %in% families) %>%
      tidyr::pivot_wider(names_from = c(family, nativeStatusCode),
                  names_prefix = "nspp_",
                  values_from = (nspp),
                  values_fill = list(nspp = 0))
  }

  # Cover by species ===========================================================
  if(!is.na(spp)){
    bysp <- full_on_cover%>%
      dplyr::group_by(site, plotID, subplotID, year) %>%
      dplyr::mutate(total_cover = sum(cover))%>%
      dplyr::ungroup() %>%
      dplyr::group_by(site, plotID, subplotID, year, scientificName) %>%
      dplyr::summarise(cover = sum(cover),
                total_cover = first(total_cover)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(rel_cover = cover/total_cover) %>%
      dplyr::ungroup() %>%
       dplyr::mutate(gen = stringr::str_split(scientificName,
                        pattern = " ",
                        simplify = TRUE)[,1],
              sp = stringr::str_split(scientificName,
                                pattern = " ",
                                simplify = TRUE)[,2],
              gen_sp = stringr::str_c(gen, " ", sp)) %>%
       dplyr::mutate(genus = stringr::str_split(scientificName,
                        pattern = " ",
                        simplify = TRUE)[,1],
              species = stringr::str_split(scientificName,
                                pattern = " ",
                                simplify = TRUE)[,2],
              gen_sp = stringr::str_c(genus, " ", species)) %>%
      dplyr::filter(gen_sp %in% spp)

    rc_sp <- bysp%>%
      dplyr::select(site, plotID, subplotID,year, gen_sp, rel_cover) %>%
      dplyr::mutate(gen_sp = stringr::str_replace(gen_sp, " ","_")) %>%
      # the following 3 lines fix some kind of duplicate problem that happened here
      # need to fix something upstream - jan 14 fixed it, need to test still
      dplyr::group_by(site, plotID, subplotID,year, gen_sp) %>%
      dplyr::summarise(rel_cover = mean(rel_cover)) %>%
      dplyr::ungroup()%>%
      tidyr::pivot_wider(names_from = gen_sp,
                  names_prefix = "rel_cover_",
                  values_from = (rel_cover),
                  values_fill = list(rel_cover = 0))

    c_sp<- bysp%>%
      dplyr::select(site, plotID, subplotID,year, gen_sp, cover) %>%
      dplyr::mutate(gen_sp = stringr::str_replace(gen_sp, " ","_")) %>%
      # the following 3 lines fix some kind of duplicate problem that happened here
      # need to fix something upstream - jan 14 fixed it, need to test still
      dplyr::group_by(site, plotID, subplotID,year, gen_sp) %>%
      dplyr::summarise(cover = mean(cover)) %>%
      dplyr::ungroup()%>%
      tidyr::pivot_wider(names_from = gen_sp,
                  names_prefix = "cover_",
                  values_from = (cover),
                  values_fill = list(cover = 0))
  }
  # by family, divided by biogeographic origin =================================
  if(!is.na(families)){
  family_stuff <- full_on_cover%>%
    dplyr::filter(nativeStatusCode %in% c("I", "N", "UNK")) %>%
    dplyr::group_by(site, plotID, subplotID, year) %>%
    dplyr::mutate(total_cover = sum(cover))%>%
    dplyr::ungroup() %>%
    dplyr::group_by(site, plotID, subplotID,year, family, nativeStatusCode) %>%
    dplyr::summarise(cover = sum(cover),
              total_cover = first(total_cover)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rel_cover = cover/total_cover) %>%
    dplyr::ungroup() %>%
    dplyr::filter(family %in% families)

  rc_ig<- family_stuff%>%
    dplyr::select(site, plotID, subplotID,year, family, nativeStatusCode,rel_cover) %>%
    dplyr::filter(nativeStatusCode == "I") %>%
    tidyr::pivot_wider(names_from = family,
                names_prefix = "rc_exotic_",
                values_from = (rel_cover),
                values_fill = list(rel_cover = 0)) %>%
    dplyr::select(-nativeStatusCode)

  rc_neg<- family_stuff%>%
    dplyr::select(site, plotID, subplotID,year, family, nativeStatusCode,rel_cover) %>%
    dplyr::filter(nativeStatusCode != "I") %>%
    tidyr::pivot_wider(names_from = family,
                       names_prefix = "rc_notexotic_",
                       values_from = (rel_cover),
                       values_fill = list(rel_cover = 0)) %>%
    dplyr::select(-nativeStatusCode)

  rc_ng<- family_stuff%>%
    dplyr::select(site, plotID, subplotID,year, family, nativeStatusCode,rel_cover) %>%
    dplyr::filter(nativeStatusCode == "N") %>%
    tidyr::pivot_wider(names_from = family,
                names_prefix = "rc_native_",
                values_from = (rel_cover),
                values_fill = list(rel_cover = 0)) %>%
    dplyr::select(-nativeStatusCode)

  c_ig<- family_stuff%>%
    dplyr::select(site, plotID, subplotID,year, family, nativeStatusCode,cover) %>%
    dplyr::filter(nativeStatusCode == "I") %>%
    tidyr::pivot_wider(names_from = family,
                names_prefix = "cover_exotic_",
                values_from = (cover),
                values_fill = list(cover = 0)) %>%
    dplyr::select(-nativeStatusCode)
  c_neg<- family_stuff%>%
    dplyr::select(site, plotID, subplotID,year, family, nativeStatusCode,cover) %>%
    dplyr::filter(nativeStatusCode != "I") %>%
    tidyr::pivot_wider(names_from = family,
                       names_prefix = "cover_notexotic_",
                       values_from = (cover),
                       values_fill = list(cover = 0)) %>%
    dplyr::select(-nativeStatusCode)
  c_ng<- family_stuff%>%
    dplyr::select(site, plotID, subplotID,year, family, nativeStatusCode,cover) %>%
    dplyr::filter(nativeStatusCode == "N") %>%
    tidyr::pivot_wider(names_from = family,
                names_prefix = "cover_native_",
                values_from = (cover),
                values_fill = list(cover = 0)) %>%
    dplyr::select(-nativeStatusCode)
  }
  # exotic  diversity and evenness ===============
  if(nrow(dplyr::filter(full_on_cover,nativeStatusCode=="I"))>0){
  vegan_friendly_div_ex <- full_on_cover %>%
    dplyr::filter(nativeStatusCode %in% c("I")) %>%
    dplyr::group_by(site, plotID, subplotID,taxonID, year) %>%
    dplyr::summarise(cover = sum(cover, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(site, plotID, subplotID, year) %>%
    tidyr::spread(taxonID, cover, fill=0) %>%
    dplyr::ungroup()

  nspp_ex <- vegan_friendly_div_ex %>%
    dplyr::select(site, plotID, subplotID,year) %>%
    dplyr::mutate(shannon_exotic = vegan::diversity(vegan_friendly_div_ex %>%
                                        dplyr::select(-site,
                                                      -plotID,
                                                      -subplotID,
                                                      -year)),
           evenness_exotic = shannon_exotic/vegan::specnumber(vegan_friendly_div_ex%>%
                                                                  dplyr::select(-site, -plotID, -subplotID, -year)),
           nspp_exotic = vegan::specnumber(vegan_friendly_div_ex %>%
                                      dplyr::select(-site,
                                                    -plotID,
                                                    -subplotID,
                                                    -year)))}else{
                                                      nspp_ex<-
                                                        template %>%
                                                        dplyr::mutate(shannon_exotic = 0,
                                                               evenness_exotic = 0,
                                                               nspp_exotic = 0)
                                                    }

  # native diversity and evenness ===========
  vegan_friendly_div_n<- full_on_cover %>%
    dplyr::filter(nativeStatusCode %in% c("N")) %>%
    dplyr::group_by(site, plotID, subplotID,taxonID, year) %>%
    dplyr::summarise(cover = sum(cover, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(taxonID = as.character(taxonID),
                  plotID = as.character(plotID)) %>%
    dplyr::group_by(site, plotID, subplotID, year) %>%
    tidyr::spread(taxonID, cover, fill=0) %>%
    dplyr::ungroup()

  nspp_n <- vegan_friendly_div_n %>%
    dplyr::select(site, plotID, subplotID,year) %>%
    dplyr::mutate(shannon_native = vegan::diversity(vegan_friendly_div_n %>%
                                               dplyr::select(-site,
                                                             -plotID,
                                                             -subplotID,
                                                             -year)),
           evenness_native = shannon_native/vegan::specnumber(vegan_friendly_div_n%>%
                                                                  dplyr::select(-site, -plotID, -subplotID, -year)),
           nspp_native = vegan::specnumber(vegan_friendly_div_n %>%
                                             dplyr::select(-site,
                                                           -plotID,
                                                           -subplotID,
                                                           -year)))


  # unknown diversity and evenness  ========================
  vegan_friendly_div_un<- full_on_cover %>%
    dplyr::filter(nativeStatusCode %in% c("UNK")) %>%
    dplyr::group_by(site, plotID, subplotID,taxonID, year) %>%
    dplyr::summarise(cover = sum(cover, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(taxonID = as.character(taxonID),
                  plotID = as.character(plotID)) %>%
    dplyr::group_by(site, plotID, subplotID, year) %>%
    tidyr::spread(taxonID, cover, fill=0) %>%
    dplyr::ungroup()

  nspp_un <- vegan_friendly_div_un %>%
    dplyr::select(site, plotID, subplotID,year) %>%
    dplyr::mutate(shannon_unknown = vegan::diversity(vegan_friendly_div_un %>%
                                               dplyr::select(-site,
                                                             -plotID,
                                                             -subplotID,
                                                             -year)),
           evenness_unknown = shannon_unknown/vegan::specnumber(vegan_friendly_div_un%>%
                                                              dplyr::select(-site, -plotID, -subplotID, -year)),
           nspp_unknown = vegan::specnumber(vegan_friendly_div_un %>%
                                             dplyr::select(-site,
                                                           -plotID,
                                                           -subplotID,
                                                           -year)))
  # not exotic diversity and evenness  ====================
  vegan_friendly_div_nex <- full_on_cover %>%
    dplyr::filter(nativeStatusCode != c("I")) %>%
    dplyr::group_by(site, plotID, subplotID,taxonID, year) %>%
    dplyr::summarise(cover = sum(cover, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(taxonID = as.character(taxonID),
                  plotID = as.character(plotID)) %>%
    dplyr::group_by(site, plotID, subplotID, year) %>%
    tidyr::spread(taxonID, cover, fill=0) %>%
    dplyr::ungroup()

  nspp_nex <- vegan_friendly_div_nex %>%
    dplyr::select(site, plotID, subplotID,year) %>%
    dplyr::mutate(shannon_notexotic = vegan::diversity(vegan_friendly_div_nex %>%
                                                dplyr::select(-site,
                                                              -plotID,
                                                              -subplotID,
                                                              -year)),
           evenness_notexotic = shannon_notexotic/vegan::specnumber(vegan_friendly_div_nex%>%
                                                                  dplyr::select(-site, -plotID, -subplotID, -year)),
           nspp_notexotic = vegan::specnumber(vegan_friendly_div_nex %>%
                                              dplyr::select(-site,
                                                            -plotID,
                                                            -subplotID,
                                                            -year)))

  # total vegan::diversity - not splitting between native status =========
  vegan_friendly_div_total <- full_on_cover %>%
    dplyr::group_by(site, plotID, subplotID, taxonID, year) %>%
    dplyr::summarise(cover = sum(cover)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(taxonID = as.character(taxonID),
           plotID = as.character(plotID)) %>%
    dplyr::filter(nchar(as.character(taxonID))>0) %>%
    dplyr::group_by(site, plotID, subplotID,year) %>%
    tidyr::spread(taxonID, cover, fill=0) %>%
    dplyr::ungroup()

  div_total <- dplyr::select(vegan_friendly_div_total, site, plotID, subplotID,year) %>%
    dplyr::mutate(shannon_total = vegan::diversity(vegan_friendly_div_total%>%
                                               dplyr::select(-site, -plotID, -subplotID, -year)),
           evenness_total = shannon_total/vegan::specnumber(vegan_friendly_div_total%>%
                                                                dplyr::select(-site, -plotID, -subplotID, -year)),
           nspp_total = vegan::specnumber(vegan_friendly_div_total%>%
                                             dplyr::select(-site, -plotID, -subplotID, -year)))

  # family diversity ===========================================================
  vegan_friendly_div_total_f <- full_on_cover %>%
    dplyr::filter(!is.na(family)) %>%
    dplyr::group_by(site, plotID, subplotID, family, year) %>%
    dplyr::summarise(cover = sum(cover)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(site, plotID, subplotID,year) %>%
    tidyr::spread(family, cover, fill=0) %>%
    dplyr::ungroup()

  div_total_f <- dplyr::select(vegan_friendly_div_total_f, site, plotID, subplotID,year) %>%
    dplyr::mutate(shannon_family = vegan::diversity(vegan_friendly_div_total_f%>%
                                              dplyr::select(-site, -plotID, -subplotID, -year)),
           evenness_family = shannon_family/vegan::specnumber(vegan_friendly_div_total_f%>%
                                               dplyr::select(-site, -plotID, -subplotID, -year)),
           nfamilies = vegan::specnumber(vegan_friendly_div_total_f%>%
                                            dplyr::select(-site, -plotID, -subplotID, -year)))

  # joining and writing out ------------------------------------------------------
  final_table <- template %>%
    dplyr::left_join(nspp_ex, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(nspp_nex, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(nspp_n, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(nspp_un, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(n_i_rel_cover, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(n_e_rel_cover, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(div_total, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(div_total_f, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::mutate(scale = scale,
                  invaded = ifelse(cover_exotic > 0, "invaded", "not_invaded"))
  if(exists("bd")){
    final_table <- final_table %>%
      dplyr::left_join(bd, by = c("site", "plotID", "subplotID", "year"))
  }
  if(!is.na(families)){
    final_table <- final_table %>%
    dplyr::left_join(rcf, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(rc_ig, by = c("site", "plotID", "subplotID", "year"))%>%
    dplyr::left_join(c_ig, by = c("site", "plotID", "subplotID", "year"))%>%
    dplyr::left_join(cf, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(rc_ng, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(c_ng, by = c("site", "plotID", "subplotID", "year"))%>%
      dplyr::left_join(rc_neg, by = c("site", "plotID", "subplotID", "year")) %>%
      dplyr::left_join(c_neg, by = c("site", "plotID", "subplotID", "year"))%>%
    dplyr::left_join(nspp_byfam, by = c("site", "plotID", "subplotID", "year"))}
  if(!is.na(spp)){
    final_table <- final_table %>%
    dplyr::left_join(rc_sp, by = c("site", "plotID", "subplotID", "year")) %>%
    dplyr::left_join(c_sp, by = c("site", "plotID", "subplotID", "year"))}

  # seems crazy, i know... but those NAs should all definitely be zero
  final_table <- final_table %>%
    dplyr::mutate_all(functions(replace(., is.na(.), 0))) %>%
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
#'
#' @examples
#'
#' # download the NEON plant diversity data for the Jornada Experimental Range
#' # raw_div<-download_plant_div(sites = "JORN")
#'
#' # convert to longform cover
#' # lf_div <- npe_longform(raw_div)
#'
#' # change all of the unknown Abutilon spp to native
#' # modified_lf_div <- change_native_status_code(lf_div, "ABUTI", "JORN", "N")
#'
#' @export
npe_change_native_status <- function(df, taxon, site, new_code){
  requireNamespace('dplyr')
  return(
    df %>%
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
#'@param dest_dir where to download the files
#'@param input to what kind of neonPlantEcology product are you appending? Can
#'be "community_matrix", "longform_cover", or "summary_info".
#'
#'
#'@export
npe_plot_centroids <- function(df,
                               dest_dir = file.path(getwd(), "tmp"),
                               type = "latlong",
                               spatial_only = T,
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

  if(input == "community_matrix") outdf <- df %>%
      tibble::rownames_to_column("plot_info") %>%
      dplyr::mutate(plotID = stringr::str_sub(plot_info,1,8)) %>%
      dplyr::left_join(neon_plots, by = "plotID")
  if(input == "longform_cover") outdf <- df %>%
      dplyr::left_join(neon_plots, by = "plotID")
  if(input == "summary_info") outdf <- df %>%
      dplyr::left_join(neon_plots, by = "plotID")

  if(spatial_only && type == "spatial") outdf <- dplyr::select(outdf, plotID )
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
#'@export
npe_plot_info <- function(comm){
  return(comm %>%
           tibble::rownames_to_column("rowname")%>%
           tidyr::separate(rowname, into = c("site", "plot", "scale", "year"), remove = F) %>%
           dplyr::mutate(plotID = stringr::str_c(site, "_", plot)) %>%
           dplyr::select(site, plot, scale, year, plotID, rowname))
}


#' get cover by species
#'
#' This extracts the cover and relative cover of all species into a data frame.
#' This is basically the same as npe_community_matrix, except the species names
#' are full, it's in a tibble format, and the plot information is not collapsed
#' into the row names.
#'
#' @export
# Cover by species ===========================================================
npe_species <- function(neon_div_object,
                        scale = "plot",
                        trace_cover = 0.5,
                        fix_unks = FALSE,
                        dissolve_years = FALSE,
                        name_cleaner = FALSE) {
  requireNamespace("tidyr")
  requireNamespace("dplyr")
  requireNamespace('vegan')
  requireNamespace("magrittr")
  # Data wrangling =============================================================

  full_on_cover <- npe_longform(neon_div_object,
                                scale = scale,
                                dissolve_years = dissolve_years,
                                fix_unks = fix_unks)
  if(name_cleaner) full_on_cover <- name_cleaner(full_on_cover)

  template <- full_on_cover %>%
    dplyr::select(site, plotID, subplotID, year)

  bysp <- full_on_cover%>%
    dplyr::group_by(site, plotID, subplotID, year) %>%
    dplyr::mutate(total_cover = sum(cover))%>%
    dplyr::ungroup() %>%
    dplyr::group_by(site, plotID, subplotID, year, scientificName) %>%
    dplyr::summarise(cover = sum(cover),
                     total_cover = first(total_cover)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rel_cover = cover/total_cover) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gen = stringr::str_split(scientificName,
                                           pattern = " ",
                                           simplify = TRUE)[,1],
                  sp = stringr::str_split(scientificName,
                                          pattern = " ",
                                          simplify = TRUE)[,2],
                  gen_sp = stringr::str_c(gen, " ", sp)) %>%
    dplyr::mutate(genus = stringr::str_split(scientificName,
                                             pattern = " ",
                                             simplify = TRUE)[,1],
                  species = stringr::str_split(scientificName,
                                               pattern = " ",
                                               simplify = TRUE)[,2],
                  gen_sp = stringr::str_c(genus, " ", species))

  rc_sp <- bysp%>%
    dplyr::select(site, plotID, subplotID,year, gen_sp, rel_cover) %>%
    dplyr::mutate(gen_sp = stringr::str_replace(gen_sp, " ","_"),
                  variable = "relative_cover") %>%
    # the following 3 lines fix some kind of duplicate problem that happened here
    # need to fix something upstream - jan 14 fixed it, need to test still
    dplyr::group_by(site, plotID, subplotID,year, gen_sp) %>%
    dplyr::summarise(rel_cover = mean(rel_cover)) %>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = gen_sp,
                       values_from = (rel_cover),
                       values_fill = list(rel_cover = 0))

  c_sp<- bysp%>%
    dplyr::select(site, plotID, subplotID,year, gen_sp, cover) %>%
    dplyr::mutate(gen_sp = stringr::str_replace(gen_sp, " ","_"),
                  variable = "absolute_cover") %>%
    # the following 3 lines fix some kind of duplicate problem that happened here
    # need to fix something upstream - jan 14 fixed it, need to test still
    dplyr::group_by(site, plotID, subplotID,year, gen_sp, variable) %>%
    dplyr::summarise(cover = mean(cover)) %>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = gen_sp,
                       values_from = (cover),
                       values_fill = list(cover = 0))

  return(list(c_sp, rc_sp))
}
