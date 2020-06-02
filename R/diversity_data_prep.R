# library(neonUtilities)
# library(tidyverse)
# library(ggpubr)
# library(vegan)
#
# options(stringsAsFactors = FALSE)
#

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
#' diversity_object <- download_plant_div(sites = "SRER")
#' @export
download_plant_div <- function(sites = "SRER"){
  require(neonUtilities)
  loadByProduct(dpID = "DP1.10058.001",
                site = sites,
                check.size = F) -> x
  return(x)
}


###############################################
# function to fix some of the unknown species #
###############################################
# source("R/unk_investigation.R")

#######################################################################
# get_longform_cover creates a long dataframe for cover from the neon #
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
#' @param neon_div_object the raw diversity data downloaded using
#' neondiversity::download_plant_div() or #' the function
#' neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".
#' @param trace_cover cover value for subplots where only occupancy was recorded
#' @param scale what level of aggregation? This can be "1m", "10m", "100m", or "plot",
#' which is the default.
#' @param fix_unks Should the unknown codes be altered with the "unk_fixer()"
#' function? Defaults to false. This requires manual investigation and editing
#' of the unk_fixer function.
#' @examples raw_div<-download_plant_div(sites = "SRER")
#' lf_div <- get_longform_cover(raw_div)
#' @export
get_longform_cover <- function(neon_div_object,
                               trace_cover=0.5,
                               scale = "plot",
                               fix_unks = FALSE){
  require(tidyverse)
  if(scale == "plot"){
    cover <- neon_div_object$div_1m2Data %>%
      mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(divDataType == "plantSpecies") %>%
      mutate(year = str_c(str_sub(endDate,1,4)))%>%
      replace_na(list(percentCover=trace_cover)) %>%
      group_by(plotID, subplotID, taxonID, year) %>%
      # dealing with the multiple bout issue by first getting the max cover
      # per sampling effort
      summarise(cover = max(percentCover),
                nativeStatusCode = first(nativeStatusCode),
                scientificName = first(scientificName),
                family = first(family)) %>%
      ungroup()  %>%
      filter(taxonID != "") %>%
      group_by(plotID, taxonID, year) %>%
      summarise(cover = sum(cover, na.rm=TRUE)/8,
                nativeStatusCode = first(nativeStatusCode),
                scientificName = first(scientificName),
                family = first(family)) %>%
      ungroup()

    traces <- neon_div_object$div_10m2Data100m2Data %>%
      mutate(endDate = as.Date(endDate)) %>%
      dplyr::filter(targetTaxaPresent == "Y") %>%
      mutate(year = str_c(str_sub(endDate,1,4)))%>%
      group_by(plotID, subplotID, taxonID, year) %>%
      summarise(cover = trace_cover,
                scientificName = first(scientificName),
                nativeStatusCode = first(nativeStatusCode),
                family = first(family)) %>%
      ungroup() %>%
      filter(taxonID != "") %>%
      group_by(plotID, taxonID, year) %>%
      summarise(cover = sum(cover, na.rm=TRUE)/12,
                nativeStatusCode = first(nativeStatusCode),
                scientificName = first(scientificName),
                family = first(family)) %>%
      ungroup()

    full_on_cover <- rbind(cover, traces) %>%
      group_by(plotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
      summarise(cover = sum(cover)) %>%
      ungroup()%>%
      mutate(site = str_sub(plotID, 1,4),
             subplotID = "plot")
    if(fix_unks) full_on_cover <- full_on_cover %>%  unk_fixer()

    return(full_on_cover)
  }

  cover8 <- neon_div_object$div_1m2Data %>%
    mutate(endDate = as.Date(endDate)) %>%
    dplyr::filter(divDataType == "plantSpecies") %>%
    mutate(year = str_c(str_sub(endDate,1,4)))%>%
    # entries in the df with no values but species was there
    # i.e. someone put the sp. code and forgot to fill in the number
    # putting as trace cover value
    replace_na(list(percentCover=trace_cover)) %>%
    mutate(endDate = as.Date(endDate)) %>%
    dplyr::filter(divDataType == "plantSpecies") %>%
    mutate(year = str_c(str_sub(endDate,1,4)))%>%
    group_by(plotID, subplotID, taxonID, year) %>%
    # dealing with the multiple bout issue by first getting the mean cover
    # per sampling effort, without aggregating, then later we'll aggregate.
    # that way, a fall-bloomer that isn't visible in spring, for example,
    # will be given its full cover value for fall, but then a species
    # that is there for both seasons will be averaged, if that makes sense
    summarise(cover = max(percentCover),
              nativeStatusCode = first(nativeStatusCode),
              scientificName = first(scientificName),
              family = first(family)) %>%
    ungroup()  %>%
    filter(taxonID != "") %>%
    mutate(subplotID = str_sub(subplotID, 1, 4))


  # 10m2,100m2 are given 0.5 (we can change later)
  # unique(x$div_10m2Data100m2Data$subplotID) # there are 12 subplots

  traces8 <- neon_div_object$div_10m2Data100m2Data %>%
    mutate(endDate = as.Date(endDate)) %>%
    dplyr::filter(targetTaxaPresent == "Y") %>%
    mutate(year = str_c(str_sub(endDate,1,4)))%>%
    group_by(plotID, subplotID, taxonID, year) %>%
    summarise(cover = trace_cover,
              scientificName = first(scientificName),
              nativeStatusCode = first(nativeStatusCode),
              family = first(family)) %>%
    ungroup() %>%
    filter(taxonID != "",
           subplotID != "31", # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
           subplotID != "32",
           subplotID != "40",
           subplotID != "41")  %>%
    mutate(subplotID = str_sub(subplotID, 1, 4))

  traces100s <- neon_div_object$div_10m2Data100m2Data %>%
    mutate(endDate = as.Date(endDate)) %>%
    dplyr::filter(targetTaxaPresent == "Y") %>%
    mutate(year = str_c(str_sub(endDate,1,4)))%>%
    group_by(plotID, subplotID, taxonID, year) %>%
    summarise(cover = trace_cover,
              scientificName = first(scientificName),
              nativeStatusCode = first(nativeStatusCode),
              family = first(family)) %>%
    ungroup() %>%
    mutate(site = str_sub(plotID, 1,4)) %>%
    filter(taxonID != "",
           subplotID == "31"| # these are the 100m2 subplots under which two 1m2 and 10m2 pairs are nested
             subplotID == "32"|
             subplotID == "40"|
             subplotID == "41")

  # aggregating at different scales ----------------------------------------------
  cover8_1m2 <- cover8 %>%
    group_by(plotID, subplotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
    summarise(cover = sum(cover)) %>%
    ungroup()%>%
    mutate(site = str_sub(plotID, 1,4))
  if(fix_unks) cover8_1m2 <- unk_fixer(cover8_1m2)

  cover8_1m2_10m2 <- rbind(cover8, traces8) %>%
    group_by(plotID,subplotID, taxonID, year, nativeStatusCode, scientificName, family) %>%
    summarise(cover = sum(cover)) %>%
    ungroup()%>%
    mutate(site = str_sub(plotID, 1,4))
  if(fix_unks) cover8_1m2_10m2<-cover8_1m2_10m2 %>%  unk_fixer()

  cover4 <- cover8_1m2_10m2 %>%
    mutate(subplotID = str_sub(subplotID, 1,2)) %>%
    rbind(traces100s) %>% # adding in the 100m2 subplots
    group_by(plotID, subplotID, year, taxonID) %>%
    summarise(cover = sum(cover), # this is summing together repeats from the rbinding
              scientificName = first(scientificName),
              nativeStatusCode = first(nativeStatusCode),
              family = first(family),
              site = first(site)) %>%
    ungroup()
  if(fix_unks) cover4 <- cover4 %>%  unk_fixer()


  if(scale == "1m") full_on_cover <- cover8_1m2
  if(scale == "10m") full_on_cover <- cover8_1m2_10m2
  if(scale == "100m") full_on_cover <- cover4

  return(full_on_cover)
}

#' Create a species abundance or occurrence matrix
#'
#' vegify creates a wide matrix of species cover or binary (presence/absence)
#' values with the plot/subplot/year as rownames. This is useful for the vegan
#' package, hence the name.
#'
#' @param neon_div_object the raw diversity data downloaded using
#' neondiversity::download_plant_div() or #' the function
#' neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".
#' @param scale what level of aggregation? This can be "1m", "10m", "100m", or "plot",
#' which is the default.
#' @param trace_cover cover value for subplots where only occupancy was recorded
#' @param fix_unks Should the unknown codes be altered with the "unk_fixer()"
#' function? Defaults to false. This requires manual investigation and editing
#' of the unk_fixer function.
#' @param binary should the matrix be converted from percent cover to binary?
#' @export
vegify <- function(neon_div_object,
                   scale="plot",
                   trace_cover = 0.5,
                   fix_unks = FALSE,
                   binary=FALSE) {
  require(tidyverse)
  if(!binary){
    return(
      neon_div_object %>%
        get_longform_cover(scale = scale, trace_cover = trace_cover, fix_unks = FALSE) %>%
        mutate(p_sp_y = paste(plotID, subplotID, year, sep = "_")) %>%
        dplyr::select(p_sp_y, taxonID, cover) %>%
        na.omit() %>% # not sure how, but there are some NA's where they shouldn't be
        pivot_wider(id_cols = p_sp_y,
                    names_from = taxonID,
                    values_from = cover,
                    values_fill = list(cover=0)) %>%
        tibble::column_to_rownames("p_sp_y") %>%
        as.data.frame()
    )
  }else{
      bin<-neon_div_object %>%
        get_longform_cover(scale = scale, trace_cover = trace_cover, fix_unks = FALSE) %>%
        mutate(p_sp_y = paste(plotID, subplotID, year, sep = "_")) %>%
        dplyr::select(p_sp_y, taxonID, cover) %>%
        na.omit() %>% # not sure how, but there are some NA's where they shouldn't be
        pivot_wider(id_cols = p_sp_y,
                    names_from = taxonID,
                    values_from = cover,
                    values_fill = list(cover=0)) %>%
        tibble::column_to_rownames("p_sp_y") %>%
        as.data.frame()

      return(ifelse(bin>0, 1,0))
  }
}


#' Get plant biodiversity information for NEON plots
#'
#' get_diversity_info calculates various biodiversity and cover indexes at the
#' plot or subplot scale at each timestep for each plot. Outputs a data frame
#' with number of species, percent cover, relative percent cover, and shannon
#' diveristy, for natives, exotics and all species. Also calculates all of these
#' metrics for the families and/or species of your choice.
#'
#' @param neon_div_object the raw diversity data downloaded using
#' neondiversity::download_plant_div() or #' the function
#' neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".
#' @param scale what level of aggregation? This can be "1m", "10m", "100m", or "plot",
#' which is the default.
#' @param trace_cover cover value for subplots where only occupancy was recorded
#' @param fix_unks Should the unknown codes be altered with the "unk_fixer()"
#' function? Defaults to false. This requires manual investigation and editing
#' of the unk_fixer function.
#' @param families Which specific families should the metrics be calculated for?
#' This can be a concatenated vector if the user want more than one family.
#' @param species Which specific species should the metrics be calculated for?
#' This can be a concatenated vector if the user want more than one species.
#' @examples
#' x <- neondiversity::download_plant_div("SRER")
#' plot_level <- get_diversity_info(neon_div_object = x, scale = "plot")
#' sp_level_1 <- get_diversity_info(x, "1m")
#' sp_level_10 <- get_diversity_info(x, "10m")
#' sp_level_100 <- get_diversity_info(x, "100m")
#' all_scales <- rbind(plot_level, sp_level_1, sp_level_10, sp_level_100)
#' @export
get_diversity_info <- function(neon_div_object,
                               scale = "plot",
                               trace_cover = 0.5,
                               fix_unks = FALSE,
                               families = "Poaceae",
                               species = "Bromus tectorum") {
  require(tidyverse)
  require(vegan)
  # Data wrangling =============================================================

  full_on_cover <- get_longform_cover(neon_div_object,
                                      scale = scale,
                                      fix_unks = fix_unks)

  # Native vs Invasive cover ===================================================

  n_i <- full_on_cover%>%
    group_by(plotID, subplotID, year) %>%
    mutate(total_cover = sum(cover))%>%
    ungroup() %>%
    group_by(plotID, subplotID,year, nativeStatusCode) %>%
    summarise(cover = sum(cover),
              total_cover = first(total_cover)) %>%
    ungroup() %>%
    mutate(rel_cover = cover/total_cover) %>%
    #filter(nativeStatusCode == "N" | nativeStatusCode == "I")%>%
    ungroup() %>%
    filter(nativeStatusCode != "")

  n_i_cover <- n_i %>%
    filter(nativeStatusCode != "" &
             nativeStatusCode != "A" &
             nativeStatusCode != "NI") %>%
    dplyr::select(plotID, subplotID,year, nativeStatusCode, cover) %>%
    pivot_wider(names_from = nativeStatusCode,
                values_from = cover,
                values_fill = list(cover = 0)) %>%
    rename(cover_native = N,
           cover_exotic = I,
           cover_unk = UNK)

  n_i_rel_cover <- n_i %>%
    filter(nativeStatusCode != ""&
             nativeStatusCode != "A" &
             nativeStatusCode != "NI") %>%
    dplyr::select(plotID, subplotID,year, nativeStatusCode, rel_cover) %>%
    pivot_wider(names_from = nativeStatusCode,
                values_from = rel_cover,
                values_fill = list(rel_cover = 0))%>%
    rename(rel_cover_native = N,
           rel_cover_exotic = I,
           rel_cover_unk = UNK)

  # Cover by family ============================================================
  byfam <- full_on_cover%>%
    group_by(plotID, subplotID, year) %>%
    mutate(total_cover = sum(cover))%>%
    ungroup() %>%
    group_by(plotID, subplotID,year, family) %>%
    summarise(cover = sum(cover),
              total_cover = first(total_cover)) %>%
    ungroup() %>%
    mutate(rel_cover = cover/total_cover) %>%
    ungroup() %>%
    filter(family %in% families)

  rcf<- byfam%>%
    dplyr::select(plotID, subplotID,year, family, rel_cover) %>%
    pivot_wider(names_from = family,
                names_prefix = "rc_",
                values_from = (rel_cover),
                values_fill = list(rel_cover = 0))

  cf<- byfam%>%
    dplyr::select(plotID, subplotID,year, family, cover) %>%
    pivot_wider(names_from = family,
                names_prefix = "cover_",
                values_from = (cover),
                values_fill = list(cover = 0))

  nspp_byfam <- full_on_cover%>%
    group_by(plotID, subplotID, year) %>%
    mutate(total_cover = sum(cover))%>%
    ungroup() %>%
    group_by(plotID, subplotID,year, family, nativeStatusCode) %>%
    summarise(nspp = length(unique(scientificName))) %>%
    ungroup() %>%
    filter(nativeStatusCode != "UNK",
           family %in% families) %>%
    pivot_wider(names_from = c(family,nativeStatusCode),
                names_prefix = "nspp_",
                values_from = (nspp),
                values_fill = list(nspp = 0))

  # Cover by species ===========================================================
  bysp <- full_on_cover%>%
    group_by(plotID, subplotID, year) %>%
    mutate(total_cover = sum(cover))%>%
    ungroup() %>%
    group_by(plotID, subplotID, year, scientificName) %>%
    summarise(cover = sum(cover),
              total_cover = first(total_cover)) %>%
    ungroup() %>%
    mutate(rel_cover = cover/total_cover) %>%
    ungroup() %>%
     mutate(gen = str_split(scientificName,
                      pattern = " ",
                      simplify = TRUE)[,1],
            sp = str_split(scientificName,
                              pattern = " ",
                              simplify = TRUE)[,2],
            gen_sp = str_c(gen, " ", sp)) %>%
     mutate(genus = str_split(scientificName,
                      pattern = " ",
                      simplify = TRUE)[,1],
            species = str_split(scientificName,
                              pattern = " ",
                              simplify = TRUE)[,2],
            gen_sp = str_c(genus, " ", species)) %>%
    filter(gen_sp %in% species)

  rc_sp <- bysp%>%
    dplyr::select(plotID, subplotID,year, gen_sp, rel_cover) %>%
    mutate(gen_sp = str_replace(gen_sp, " ","_")) %>%
    # the following 3 lines fix some kind of duplicate problem that happened here
    # need to fix something upstream - jan 14 fixed it, need to test still
    group_by(plotID, subplotID,year, gen_sp) %>%
    summarise(rel_cover = mean(rel_cover)) %>%
    ungroup()%>%
    pivot_wider(names_from = gen_sp,
                names_prefix = "rc_",
                values_from = (rel_cover),
                values_fill = list(rel_cover = 0))

  c_sp<- bysp%>%
    dplyr::select(plotID, subplotID,year, gen_sp, cover) %>%
    mutate(gen_sp = str_replace(gen_sp, " ","_")) %>%
    # the following 3 lines fix some kind of duplicate problem that happened here
    # need to fix something upstream - jan 14 fixed it, need to test still
    group_by(plotID, subplotID,year, gen_sp) %>%
    summarise(cover = mean(cover)) %>%
    ungroup()%>%
    pivot_wider(names_from = gen_sp,
                names_prefix = "cover_",
                values_from = (cover),
                values_fill = list(cover = 0))

  # by family, divided by biogeographic origin =================================
  exotic_grass <- full_on_cover%>%
    group_by(plotID, subplotID, year) %>%
    mutate(total_cover = sum(cover))%>%
    ungroup() %>%
    group_by(plotID, subplotID,year, family, nativeStatusCode) %>%
    summarise(cover = sum(cover),
              total_cover = first(total_cover)) %>%
    ungroup() %>%
    mutate(rel_cover = cover/total_cover) %>%
    ungroup() %>%
    filter(family %in% families)

  rc_ig<- exotic_grass%>%
    dplyr::select(plotID, subplotID,year, family, nativeStatusCode,rel_cover) %>%
    filter(nativeStatusCode == "I") %>%
    pivot_wider(names_from = family,
                names_prefix = "rc_exotic_",
                values_from = (rel_cover),
                values_fill = list(rel_cover = 0)) %>%
    dplyr::select(-nativeStatusCode)

  rc_ng<- exotic_grass%>%
    dplyr::select(plotID, subplotID,year, family, nativeStatusCode,rel_cover) %>%
    filter(nativeStatusCode == "N") %>%
    pivot_wider(names_from = family,
                names_prefix = "rc_native_",
                values_from = (rel_cover),
                values_fill = list(rel_cover = 0)) %>%
    dplyr::select(-nativeStatusCode)

  c_ig<- exotic_grass%>%
    dplyr::select(plotID, subplotID,year, family, nativeStatusCode,cover) %>%
    filter(nativeStatusCode == "I") %>%
    pivot_wider(names_from = family,
                names_prefix = "cover_exotic_",
                values_from = (cover),
                values_fill = list(cover = 0)) %>%
    dplyr::select(-nativeStatusCode)

  c_ng<- exotic_grass%>%
    dplyr::select(plotID, subplotID,year, family, nativeStatusCode,cover) %>%
    filter(nativeStatusCode == "N") %>%
    pivot_wider(names_from = family,
                names_prefix = "cover_native_",
                values_from = (cover),
                values_fill = list(cover = 0)) %>%
    dplyr::select(-nativeStatusCode)

  # diversity indexes splitting between native status ==========================
  vegan_friendly_div <- full_on_cover %>%
    group_by(plotID, subplotID,taxonID, year, nativeStatusCode) %>%
    summarise(cover = sum(cover, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(taxonID = as.character(taxonID),
           plotID = as.character(plotID),
           nativeStatusCode = as.character(nativeStatusCode)) %>%
    filter(nchar(as.character(taxonID))>0,
           nativeStatusCode != "",
           nativeStatusCode != "A",
           nativeStatusCode != "NI") %>%
    group_by(plotID, subplotID, year, nativeStatusCode) %>%
    spread(taxonID, cover, fill=0) %>%
    ungroup()

  # note to self - hard code! gotta fix it
  vegan_friendly_div$shannon = diversity(vegan_friendly_div[5:ncol(vegan_friendly_div)])
  vegan_friendly_div$nspp = specnumber(vegan_friendly_div[5:ncol(vegan_friendly_div)])

  nspp <- vegan_friendly_div %>%
    dplyr::select(plotID, subplotID,year, nativeStatusCode, nspp) %>%
    pivot_wider(names_from = nativeStatusCode,
                values_from = nspp,
                values_fill = list(nspp=0)) %>%
    rename(nspp_native = N, nspp_exotic=I, nspp_unk = UNK)

  shannon <- vegan_friendly_div %>%
    dplyr::select(plotID, subplotID,year, nativeStatusCode, shannon) %>%
    pivot_wider(names_from = nativeStatusCode,
                values_from = shannon,
                values_fill = list(shannon=0)) %>%
    rename(shannon_native = N, shannon_exotic=I, shannon_unk = UNK)

  # total diversity - not splitting between native status
  vegan_friendly_div_total <- full_on_cover %>%
    group_by(plotID, subplotID, taxonID, year) %>%
    summarise(cover = sum(cover)) %>%
    ungroup() %>%
    mutate(taxonID = as.character(taxonID),
           plotID = as.character(plotID)) %>%
    filter(nchar(as.character(taxonID))>0) %>%
    group_by(plotID, subplotID,year) %>%
    spread(taxonID, cover, fill=0) %>%
    ungroup()

  # note to self - fix that hard-coding!
  div_total <- dplyr::select(vegan_friendly_div_total, plotID, subplotID,year)
  div_total$shannon_total = diversity(vegan_friendly_div_total[4:ncol(vegan_friendly_div_total)])
  div_total$nspp_total = specnumber(vegan_friendly_div_total[4:ncol(vegan_friendly_div_total)])

  # joining and writing out ------------------------------------------------------
  final_table <- left_join(nspp, shannon, by = c("plotID", "subplotID", "year")) %>%
    left_join(n_i_cover, by = c("plotID", "subplotID","year")) %>%
    left_join(n_i_rel_cover, by = c("plotID", "subplotID", "year")) %>%
    left_join(div_total, by = c("plotID", "subplotID", "year"))%>%
    left_join(rcf, by = c("plotID", "subplotID", "year")) %>%
    left_join(rc_ig, by = c("plotID", "subplotID", "year"))%>%
    left_join(c_ig, by = c("plotID", "subplotID", "year"))%>%
    left_join(cf, by = c("plotID", "subplotID", "year")) %>%
    left_join(rc_ng, by = c("plotID", "subplotID", "year")) %>%
    left_join(c_ng, by = c("plotID", "subplotID", "year")) %>%
    left_join(rc_sp, by = c("plotID", "subplotID", "year")) %>%
    left_join(c_sp, by = c("plotID", "subplotID", "year")) %>%
    left_join(nspp_byfam, by = c("plotID", "subplotID", "year")) %>%
    mutate(site = str_sub(plotID, 1,4),
           scale = scale,
           invaded = if_else(cover_exotic > 0, "invaded", "not_invaded"))%>%
    mutate(scale = factor(scale, levels = c("1m","10m","100m", "plot")))

  # seems crazy, i know... but those NAs should all definitely be zero
  final_table[is.na(final_table)] <- 0

  return(final_table)
}

# using the functions ==========================================================
#
# plot_level <- get_diversity_info(neon_div_object = x, scale = "plot")
# sp_level_1 <- get_diversity_info(x, "1m")
# sp_level_10 <- get_diversity_info(x, "10m")
# sp_level_100 <- get_diversity_info(x, "100m")
# all_scales <- rbind(plot_level, sp_level_1, sp_level_10, sp_level_100)


# write_csv(plot_level, "data/plot_level_diversity.csv")
# write_csv(sp_level_1, "data/subplot_level_diversity_1.csv")
# write_csv(sp_level_10, "data/subplot_level_diversity_10.csv")
# write_csv(sp_level_100, "data/subplot_level_diversity_100.csv")

# # for example
# test <- get_diversity_info(x,"plot",
#                            families = c("Poaceae", "Brassicaceae"),
#                            species = c("Bromus tectorum", "Poa secunda"))
# summary(test)

