# reclassifying unknowns. meant to by sourced from diverisity_data_prep.R
# There is a flora specifically for the Jornada
# there is probably a plant list for the santa rita experimental range on seinet
# but that site is down at the moment

# # side tangent - to look at remaining unks do something like this
# unks <- full_on_cover %>%
#   filter(nativeStatusCode == "UNK") %>%
#   select(taxonID, plotID, family, scientificName) %>%
#   mutate(site = str_sub(plotID, 1,4)) %>%
#   group_by(site, taxonID) %>%
#   summarise(family = first(family),
#             scientificName = paste(unique(scientificName))) %>%
#   ungroup()


#' Manually change unknown native status codes
#'
#' Usually at any NEON site, there will be plenty of unknown plants for various
#' reasons. Often, it is possible to use the information that is available to
#' at least determine if a particular unknown species is native or introduced.
#' For example, it might be documented in the local flora that all plants of a
#' certain family or genus are native. This function allows the user to manually
#' change those codes. However it must be manually edited at this point.
#'
#' @export
unk_fixer <- function(df){
  return(
    df %>%
      mutate(nativeStatusCode = replace(nativeStatusCode, taxonID == "ABUTI", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ASTRA", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CORYP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CRYPT" &
                                 site == "JORN", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CRYPTSPP" &
                                 site == "JORN", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CUSCU", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ERIOG", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "EUPHO" &
                                 site == "JORN", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "IPOMO", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "LINUM" &
                                          site == "JORN", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "MALVAC" &
                                          site == "JORN", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "OPUNT", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "PANIC" &
                                          site == "JORN", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "PROSO", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "SENNA", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "SPORO", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "SPOROSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ARTEM", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ASTRASPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ATRIPSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CACTAC", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CACTACSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "DRABA", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "DRABASPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ERIOGSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ELYMU", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "OENOT", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "OENOTSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "OPUNTSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "PEDIO", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "PEDIOSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "PHYSA2", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "PHYSASPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "POLEMO", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "POLEMOSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "SCLER10", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "TOWNS", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "TOWNSSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ACHNA", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ALLIU", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ARENA", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ATRIP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "BROMU", "I"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CASTI2", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CASTI2SPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CREPI", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "CREPISPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "EPILO", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ERIGE2", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ERIGE2SPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "OROBA", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "OROBASPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "SCROPH", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "SCROPHSPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "TRAGOSPP", "I"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "AMBRO", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "AYENI", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "BIDEN", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "ECHIN3", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "LYCIU", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "SIDA", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "SIDASPP", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "SISYM", "N"),
             nativeStatusCode = replace(nativeStatusCode, taxonID == "THYMO", "N")
             )
  )
# full_on_cover[full_on_cover$taxonID == "IPOMOSPP","nativeStatusCode"] <-"N"
# full_on_cover[full_on_cover$taxonID == "LYCIU","nativeStatusCode"] <-"N"        # Lycium
# full_on_cover[full_on_cover$taxonID == "SIDA","nativeStatusCode"] <-"N"
# full_on_cover[full_on_cover$taxonID == "SIDASPP","nativeStatusCode"] <-"N"
# full_on_cover[full_on_cover$taxonID == "SISYM","nativeStatusCode"] <-"I"
# full_on_cover[full_on_cover$taxonID == "THYMO","nativeStatusCode"] <-"N"
}




# # investigating erodium unk
# full_on_cover[full_on_cover$taxonID == "ERODI",]
# full_on_cover[full_on_cover$plotID == "SRER_013",] %>%
#   filter(str_sub(taxonID,1,2) == "ER" & family == "Geraniaceae")
# # need to check the SRER plant list -- these are usually not hard to figure out
# # looking at the unknown lepidiums in srer
# x$div_1m2Data[str_sub(x$div_1m2Data$scientificName,1,4) == "Lepi",] %>%
#   filter(str_sub(plotID, 1,4) == "SRER") %>%
#   dplyr::select(scientificName) %>%
#   unique()
# # so, there's densiflorum and lasiocarpum to choose from. lasiocarpum has hair,
# # densiflorum can be glabrous (according to FNA)
