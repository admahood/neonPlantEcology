[![DOI](https://zenodo.org/badge/268667267.svg)](https://zenodo.org/badge/latestdoi/268667267)

# neondiversity

This package contains scripts for processing plant species cover and occurrence data from NEON into formats that are easy to use, particularily with the `vegan` package. The main functions are `vegify` which creates vegan-friendly species occurrence and abundance matrices, and `get_diversity_info`, which creates summary statistics by plot (or subplot) on the diversity, cover, relative cover, and number of species of natives, non-natives, and members of the family or species of your choosing.

## Installation Instructions

`devtools::install_github("admahood/neondiversity")`

## Use

### Download

First, download some data using the `download_plant_div` function. It defaults to download the diversity data for the Santa Rita Experimental Range in Domain 14 in Arizona. The `sites` arguement in this function is used to specify which site (or sites) you wish to download. A list of field sites can be found [here](https://www.neonscience.org/field-sites/field-sites-map/list).

`sites_adam_worked_at_back_in_the_day <- download_plant_div(sites = c("SRER", "JORN")`

The output is a list of four things. The first two are of most interest. The first list item is the abundances observed in the 1m^2 subplots. The second list item is the occurrences observed for the 10m^2 and 100m^2 subplots.

### get_longform_cover

This is mainly a helper function for `vegify` and `get_diversity_info`. But if all you want is a longform dataframe of the percent cover of each species in each plot or subplot, this is the function for you.

### vegify

This function converts the diversity object downloaded from NEON into a matrix of either abundances (percent cover from 0-100) or occurrences (0 or 1), at the scale of your choosing (1m^2, 10m^2, 100m^2, or 400m^2, which is a whole plot).

`species_occurrence_matrix <- vegify(sites_adam_worked_at_back_in_the_day, binary=TRUE)`

### get_diversity_info

This function supplies the user with summary statistics and cover data for exotic plants, native plants, and all plants. 


|Variables             |
|:---------------------|
|plotID                |
|subplotID             |
|year                  |
|nspp_native           |
|nspp_exotic           |
|nspp_unk              |
|shannon_native        |
|shannon_exotic        |
|shannon_unk           |
|cover_native          |
|cover_exotic          |
|cover_unk             |
|rel_cover_native      |
|rel_cover_exotic      |
|rel_cover_unk         |
|shannon_total         |
|nspp_total            |
|rc_Poaceae            |
|rc_exotic_Poaceae     |
|cover_exotic_Poaceae  |
|cover_Poaceae         |
|rc_native_Poaceae     |
|cover_native_Poaceae  |
|rc_Bromus_tectorum    |
|cover_Bromus_tectorum |
|nspp_Poaceae_N        |
|nspp_Poaceae_I        |
|site                  |
|scale                 |
|invaded               |


|   |   plotID        | subplotID       |    year         | nspp_native  | nspp_exotic   |   nspp_unk    |shannon_native |shannon_exotic | shannon_unk   | cover_native  | cover_exotic   |  cover_unk     |rel_cover_native |rel_cover_exotic |rel_cover_unk   |shannon_total |  nspp_total  |  rc_Poaceae   |rc_exotic_Poaceae |cover_exotic_Poaceae |cover_Poaceae   |rc_native_Poaceae |cover_native_Poaceae |rc_Bromus_tectorum |cover_Bromus_tectorum |nspp_Poaceae_N |nspp_Poaceae_I |    site         | scale    |  invaded        |
|:--|:----------------|:----------------|:----------------|:-------------|:--------------|:--------------|:--------------|:--------------|:--------------|:--------------|:---------------|:---------------|:----------------|:----------------|:---------------|:-------------|:-------------|:--------------|:-----------------|:--------------------|:---------------|:-----------------|:--------------------|:------------------|:---------------------|:--------------|:--------------|:----------------|:---------|:----------------|
|   |Length:11625     |Length:11625     |Length:11625     |Min.   : 0.00 |Min.   : 0.000 |Min.   : 0.000 |Min.   :0.0000 |Min.   :0.0000 |Min.   :0.0000 |Min.   :  0.00 |Min.   :  0.000 |Min.   : 0.0000 |Min.   :0.0000   |Min.   :0.00000  |Min.   :0.00000 |Min.   :0.000 |Min.   : 1.00 |Min.   :0.0000 |Min.   :0.00000   |Min.   :  0.000      |Min.   :  0.000 |Min.   :0.00000   |Min.   :  0.000      |Min.   :0.00000    |Min.   : 0.0000       |Min.   : 0.000 |Min.   :0.000  |Length:11625     |1m  :4326 |Length:11625     |
|   |Class :character |Class :character |Class :character |1st Qu.: 4.00 |1st Qu.: 0.000 |1st Qu.: 0.000 |1st Qu.:0.6931 |1st Qu.:0.0000 |1st Qu.:0.0000 |1st Qu.:  3.50 |1st Qu.:  0.000 |1st Qu.: 0.0000 |1st Qu.:0.7000   |1st Qu.:0.00000  |1st Qu.:0.00000 |1st Qu.:1.099 |1st Qu.: 5.00 |1st Qu.:0.1282 |1st Qu.:0.00000   |1st Qu.:  0.000      |1st Qu.:  1.000 |1st Qu.:0.05571   |1st Qu.:  0.500      |1st Qu.:0.00000    |1st Qu.: 0.0000       |1st Qu.: 1.000 |1st Qu.:0.000  |Class :character |10m :4380 |Class :character |
|   |Mode  :character |Mode  :character |Mode  :character |Median : 7.00 |Median : 1.000 |Median : 1.000 |Median :1.3750 |Median :0.0000 |Median :0.0000 |Median :  8.00 |Median :  0.500 |Median : 0.5000 |Median :0.8750   |Median :0.02256  |Median :0.03448 |Median :1.644 |Median : 8.00 |Median :0.3000 |Median :0.00000   |Median :  0.000      |Median :  2.500 |Median :0.18182   |Median :  1.500      |Median :0.00000    |Median : 0.0000       |Median : 2.000 |Median :0.000  |Mode  :character |100m:2333 |Mode  :character |
|   |NA               |NA               |NA               |Mean   : 9.76 |Mean   : 1.649 |Mean   : 1.696 |Mean   :1.3518 |Mean   :0.3322 |Mean   :0.3209 |Mean   : 16.17 |Mean   :  1.678 |Mean   : 0.9561 |Mean   :0.7893   |Mean   :0.12290  |Mean   :0.08743 |Mean   :1.641 |Mean   :11.51 |Mean   :0.3485 |Mean   :0.06302   |Mean   :  1.123      |Mean   :  5.632 |Mean   :0.26207   |Mean   :  4.352      |Mean   :0.03536    |Mean   : 0.3805       |Mean   : 2.765 |Mean   :0.465  |NA               |plot: 586 |NA               |
|   |NA               |NA               |NA               |3rd Qu.:12.00 |3rd Qu.: 3.000 |3rd Qu.: 3.000 |3rd Qu.:1.9342 |3rd Qu.:0.6931 |3rd Qu.:0.6931 |3rd Qu.: 20.00 |3rd Qu.:  1.500 |3rd Qu.: 1.0000 |3rd Qu.:0.9688   |3rd Qu.:0.15385  |3rd Qu.:0.11538 |3rd Qu.:2.164 |3rd Qu.:14.00 |3rd Qu.:0.5227 |3rd Qu.:0.05714   |3rd Qu.:  0.500      |3rd Qu.:  6.417 |3rd Qu.:0.41226   |3rd Qu.:  5.000      |3rd Qu.:0.03030    |3rd Qu.: 0.5000       |3rd Qu.: 4.000 |3rd Qu.:1.000  |NA               |NA        |NA               |
|   |NA               |NA               |NA               |Max.   :86.00 |Max.   :13.000 |Max.   :18.000 |Max.   :4.1494 |Max.   :2.0621 |Max.   :2.5239 |Max.   :200.00 |Max.   :100.500 |Max.   :42.0000 |Max.   :1.0000   |Max.   :1.00000  |Max.   :1.00000 |Max.   :4.175 |Max.   :98.00 |Max.   :1.0000 |Max.   :1.00000   |Max.   :100.500      |Max.   :126.000 |Max.   :1.00000   |Max.   :126.000      |Max.   :1.00000    |Max.   :85.0000       |Max.   :19.000 |Max.   :4.000  |NA               |NA        |NA               |

### unk_fixer
