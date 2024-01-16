[![DOI](https://zenodo.org/badge/268667267.svg)](https://zenodo.org/badge/latestdoi/268667267)

# neonPlantEcology

The National Ecological Observatory Network (NEON) collects long-term ecological monitoring data on myriad ecosystem components, including plant diversity. Plant diversity data is collected at yearly or sub-yearly time steps, depending on the ecology of the site, at plots that are constructed with a nested design. The spatial and temporal nature of the sampling design, and the resulting storage of the data, may not be straightforward to the average end user. Thus, I created an R package for transforming this raw data product into forms that are easy to use for ecologists. The following is a short explanation of the data transformation process.

This package contains scripts for processing plant species cover and occurrence data from NEON into formats that are easy to use, particularily with the `vegan` package. The main functions are `npe_community_matrix` which creates vegan-friendly species occurrence and abundance matrices, and `npe_summary`, which creates summary statistics by site, plot or subplot on the diversity, cover, relative cover, and number of species of natives, non-natives, and members of the family or species of your choosing.

## Installation Instructions

#### development version

`devtools::install_github("admahood/neonPlantEcology")`

#### cran version

`install.packages("neonPlantEcology")`

## Use

### npe_download

First, download some data using the `download_plant_div` function. It defaults to download the diversity data for the Santa Rita Experimental Range in Domain 14 in Arizona. The `sites` arguement in this function is used to specify which site (or sites) you wish to download. A list of field sites can be found [here](https://www.neonscience.org/field-sites/explore-field-sites).

`sites <- npe_download_plant_div(sites = c("SRER", "JORN")`

The output is a list of four things. The first two are of most interest. The first list item is the abundances observed in the 1m<sup>2</sup> subplots. The second list item is the occurrences observed for the 10m<sup>2</sup> and 100m<sup>2</sup> subplots.

The function `npe_download` is more generalizaed, and can be used to download any data product from neon.

### npe_community_matrix

This function converts the diversity object downloaded from NEON into a matrix of either abundances (percent cover from 0-100) or occurrences (0 or 1), at the scale of your choosing (1m<sup>2</sup>, 10m<sup>2</sup>, 100m<sup>2</sup>, or 400m<sup>2</sup>, which is a whole plot). It can also take as an input the output from the `divStack` function in the neonPlants package.

`species_occurrence_matrix <- npe_community_matrix(sites, binary=TRUE)`

### npe_summary (NOTE: formerly npe_diversity_info)

`npe_summary` calculates various biodiversity and cover indexes at the plot or subplot scale for each year for each plot. Outputs a data frame with number of species, percent cover, relative percent cover, and shannon diveristy, for natives, exotics and all species. Also calculates all of these metrics for the families and/or species of your choice.

| Variable | Description | Additional arguments |
|:---------|-------------|-------------|
|shannon_<exotic/native/unknown/total> | Shannon-Weaver diversity of exotic/native/unknown or all species| Given by default |
|evennness_<exotic/native/unknown/total> | Pielou's evenness | |
|nspp_<exotic/native/unknown/total> | number of species||
|cover_<exotic/native/unknown/total> | Absolute cover as measured by technicians||
|rel_cover_<exotic/native/unknown/total> | Relative cover - the absoulte cover divided by the total cover of all species||
|nfamilies| number of families ||
|shannon_family| Shannon-Weaver diversity, but aggregated by family instead of species||
|evenness_family| Pielou's evenness, but aggregated by family instead of species||
|scale | the scale of aggregation (1m, 10m, 100m, plot or site)||
|invaded | is there at least one exotic species present?||
|turnover | species turnover according to vegan::nestedbetajac() (Baselga 2012) |betadiversity = TRUE, scale = c("plot", "site")|
|nestedness | nestedness according to vegan::nestedbetajac() (Baselga 2012) |betadiversity = TRUE, scale = c("plot", "site")|

### npe_longform

This is really the meat of the package. It is used as a helper function for `npe_community_matrix` and `npe_summary`. In many cases the end user will not need to deal with it, but if all you want is a longform dataframe of the percent cover of each species in each plot or subplot, this is the function for you. 

After the data is downloaded from the site(s) of interest, the user needs to make three decisions about how the data is processed. First, the scale of analysis. The NEON plots have a nested design, in which each 400m2 plot is divided into 4 100m2 subplots, each of which has two nested 10m2 and 1m2 subplots (Barnett et al. 2019). Cover is estimated in the 1m2 subplots, occurrence for those species which are not present in the 1m2 subplot is recorded in each 10m2 subplot, and finally occurrence is then recorded in the 100m2 subplots for those species that do not occur in the 1m2 and 10m2 subplots that are nested within. The user can decide whether they want to analyze only the 1m2 subplots, the 1 and 10m2 subplots, 1, 10 and 100m2 or the whole plot scale. Next, the user decides what the trace cover estimate would be for the occurrence data. This is set at a default value of 0.5 percent. Lastly, there is an argument whether to fix the native status codes for unknown species.

There are a few decisions I made that are not provided as arguments to the functions. First, some sites are sampled in multiple bouts, due to the temporal heterogeneity in the life cycles of different plants. For example, in the Jornada Experimental Range, (site abbreviation “JORN”) there spring and fall blooms, and so plants are sampled around peak biomass for each of those seasons. I reasoned that this would lead to many herbaceous species that would be undetected in one season and at full biomass in the other. Thus, in order to capture a cover estimate that corresponds to peak biomass, we take the maximum cover estimate between sampling bouts. 

In some cases, the value for the cover of a species in a 1m2 plot was NA. I assumed these were cases where the field botanist recorded the species present and forgot to estimate the cover for that species, and returned to the office with a blank space on the data sheet where the cover should be. When cases like this were encountered, we replaced those NA values with the trace cover value.

These decisions are given as arguments to the function call. The function aggregates the cover to the appropriate scale. At the whole plot scale, for example, first the 1m2 subplots are aggregated by taking the intrayear max per subplot, then taking the sum of those values for the entire plot divided by the number of subplots. The 10 and 100m2 subplots are all given the trace cover value, summed, and divided by the number of subplots. Then the cover data frame is combined with the data frame of trace cover values, and summed by plot.

### npe_change_native_status

There are a many instances where plants are classified as unknown, but still have a certain level of known-ness. For example, an unknown Opuntia sp. may have the genus and family recorded, but still be listed as unknown for the nativeStatusCode. It is possible to look at local flora and determine that while the exact species is unknown, the plant is very likely to be native or non-native. For example, a site may have 5 species of Opuntia that are all native. In this case, if the user cares about the native status, we provide a script that can be edited manually (unk_investigation.R) to address this concern. I did a very basic first try for the sites we were interested in examining while I was writing these functions. All of the above functions have an arguement, `fix_unks`, that defaults to `FALSE`. If it is set to `TRUE`, the `unk_fixer` function is called from the unk_investigation.R script and native status codes are changed according to the script. This is still very much in the beta stage. Another thought on dealing with unknown species is using a model that incorporates that uncertainty developed recently by Anna I Spiers (2021).

# References

Barnett, D. T., Adler, P. B., Chemel, B. R., Duffy, P. A., Enquist, B. J., Grace, J. B., … Vellend, M. (2019). The plant diversity sampling design for The National Ecological Observatory Network. Ecosphere, 10(2), e02603. https://doi.org/10.1002/ecs2.2603

Baselga, A. (2012). The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. Global Ecol. Biogeogr. 21, 1223–1232.

Spiers, Anna I. and Royle, J. Andrew and Torrens, Christa L. and Joseph, Maxwell B. 2021. Estimating occupancy dynamics and encounter rates with species misclassification: a semi-supervised individual-level approach. Biorxiv, http://biorxiv.org/lookup/doi/10.1101/2021.03.17.433917


