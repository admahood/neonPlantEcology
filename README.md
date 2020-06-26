[![DOI](https://zenodo.org/badge/268667267.svg)](https://zenodo.org/badge/latestdoi/268667267)

# neondiversity

This package contains scripts for processing plant species cover and occurrence data from NEON into formats that are easy to use, particularily with the `vegan` package. The main functions are `vegify` which creates vegan-friendly species occurrence and abundance matrices, and `get_diversity_info`, which creates summary statistics by plot (or subplot) on the diversity, cover, relative cover, and number of species of natives, non-natives, and members of the family or species of your choosing.

## Installation Instructions

`devtools::install_github("admahood/neondiversity")`

## Use

### Download

First, download some data using the `download_plant_div` function. It defaults to download the diversity data for the Santa Rita Experimental Range in Domain 14 in Arizona. The `sites` arguement in this function is used to specify which site (or sites) you wish to download. A list of field sites can be found [here](https://www.neonscience.org/field-sites/field-sites-map/list).

`sites_adam_worked_at_back_in_the_day <- download_plant_div(sites = c("SRER", "JORN")`

The output is a list of four things. The first two are of most interest. The first list item is the abundances observed in the 1m<sup>2</sup> subplots. The second list item is the occurrences observed for the 10m<sup>2</sup> and 100m<sup>2</sup> subplots.

### get_longform_cover

This is mainly a helper function for `vegify` and `get_diversity_info`. But if all you want is a longform dataframe of the percent cover of each species in each plot or subplot, this is the function for you.

### vegify

This function converts the diversity object downloaded from NEON into a matrix of either abundances (percent cover from 0-100) or occurrences (0 or 1), at the scale of your choosing (1m<sup>2</sup>, 10m<sup>2</sup>, 100m<sup>2</sup>, or 400m<sup>2</sup>, which is a whole plot).

`species_occurrence_matrix <- vegify(sites_adam_worked_at_back_in_the_day, binary=TRUE)`

### get_diversity_info

`get_diversity_info` calculates various biodiversity and cover indexes at the plot or subplot scale for each year for each plot. Outputs a data frame with number of species, percent cover, relative percent cover, and shannon diveristy, for natives, exotics and all species. Also calculates all of these metrics for the families and/or species of your choice.

### unk_fixer

There are a many instances where plants are classified as unknown, but still have a certain level of known-ness. For example, an unknown Opuntia sp. may have the genus and family recorded, but still be listed as unknown for the nativeStatusCode. It is possible to look at local flora and determine that while the exact species is unknown, the plant is very likely to be native or non-native. For example, a site may have 5 species of Opuntia that are all native. In this case, if the user cares about the native status, we provide a script that can be edited manually (unk_investigation.R) to address this concern. I did a very basic first try for the sites we were interested in examining while I was writing these functions. All of the above functions have an arguement, `fix_unks`, that defaults to `FALSE`. If it is set to `TRUE`, the `unk_fixer` function is called from the unk_investigation.R script and native status codes are changed according to the script. This is still very much in the beta stage.
