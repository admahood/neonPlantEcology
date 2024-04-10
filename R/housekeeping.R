# package-wide setup
packageStartupMessage("Welcome to neonPlantEcology!\n
                      Type vignette('neonPlantEcology') for a quick orientation.\n
                      Please report any bugs, other issues, or desires for future capabilities to https://github.com/admahood/neonPlantEcology/issues")


# solving the note on unused imports (we do in fact use these packages for)
ignore_unused_imports <- function(){
  utils::globalVariables
  ggplot2::ggplot
  ggpubr::ggarrange
  sf::st_read
}

# solving the note on global variables (a symptom of using dplyr a lot)
utils::globalVariables(c("nativeStatusCode", "taxonID", "rowname", "site", "lut",
                         "eventID", "plotID", "subplotID", "targetTaxaPresent",
                         "present", "na.omit", "p_sp_y", "cover", "total_cover",
                         "rel_cover", "contains", "family", "scientificName",
                         "nspp", "shannon_exotic", "shannon_native", "shannon_total",
                         "shannon_notexotic", "shannon_family", "cover_exotic",
                         "endDate", "divDataType", "otherVariables", "percentCover",
                         "heightPlantSpecies", "all_na", "height", "str_detect",
                         "download.file", "unzip", "plot_info", "latitude", "longitude",
                         "domainName", "domainNumb", "siteType", "ai_class", "koppen_coarse",
                         "koppen_fine", "shannon_unknown", "siteID", "data",
                         "plot_centroids", 'bout'))
