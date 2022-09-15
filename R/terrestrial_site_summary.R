# getting summary statistics

library(tidyverse)
library(neondiveRsity)
library(neonUtilities)

# getting all the sites
sites <- read_csv("data/NEON_Field_Site_Metadata_20220412.csv")

site_codes <- sites %>%
  filter(field_site_type %in% c("Gradient Terrestrial", "Core Terrestrial")) %>%
  pull(field_site_id)

# downloading all data
if(!file.exists("data/all_div_site.Rda")){div<- list()
  for(i in 1:length(site_codes)){
    print(site_codes[i])
    div[[i]] <- neondiveRsity::download_plant_div(sites = site_codes[i]) %>%
      neondiveRsity::get_diversity_info(scale = "site")
  }

all <- bind_rows(div)
save(all, file= "data/all_div_site.Rda")
}else{load("data/all_div_site.Rda")}

# interpreting all the data

unaggregated<- all %>%
  dplyr::select(site, year, nspp_exotic, rel_cover_exotic,
                nspp_unknown, rel_cover_unknown) %>%
  pivot_longer(cols = c(nspp_exotic, rel_cover_exotic,
                        nspp_unknown, rel_cover_unknown)) %>%
  ggplot(aes(x= value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free") +
  theme_classic()
ggsave(unaggregated, filename = "summary_all_terrestrial_sites.png")

caterpillar_data <-
  all %>%
  group_by(site) %>%
  mutate(year = as.numeric(year)) %>%
  summarise(min_year = min(year),
            max_year = max(year),
            max_nspp_exotic = max(nspp_exotic),
            min_nspp_exotic = min(nspp_exotic),
            max_nspp_unknown = max(nspp_unknown),
            min_nspp_unknown = min(nspp_unknown),
            max_rel_cover_exotic = max(rel_cover_exotic),
            min_rel_cover_exotic = min(rel_cover_exotic),
            max_rel_cover_unknown = max(rel_cover_unknown),
            min_rel_cover_unknown = min(rel_cover_unknown)) %>%
  ungroup() %>%
  pivot_longer(cols = names(.)[2:ncol(.)]) %>%
  mutate(stat = str_sub(name, 1,3),
         property = str_remove_all(name, c("max_", "min_")),
         property = ifelse(property == "min_year", "year", property),
         property = ifelse(property == "max_year", "year", property)) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = c(stat), values_from = value)

order_d<-caterpillar_data %>%
  group_by(site) %>%
  filter(property == "rel_cover_exotic") %>%
  dplyr::select(site, orderval=max)

caterpillar <-
  caterpillar_data %>%
  filter(property != "year") %>%
  left_join(order_d) %>%
  ggplot(aes(y = reorder(site, orderval))) +
  geom_segment(aes(x = min, xend=max, yend = site), lwd =1) +
  facet_grid(~property, scales = "free") +
  theme_classic() +
  theme(panel.background = element_rect(fill=NA, color="black"),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_line(color = "grey", linetype = 2)) +
  xlab("Range of Values") +
  ggtitle("NEON Terrestrial Sites");caterpillar
# dir.create("figs")
ggsave(caterpillar, filename = "figs/caterpillar_plot.png", height=7,width = 10)


all %>%
  group_by(site) %>%
  summarise(nspp_exotic = max(nspp_exotic),
            nspp_unknown = max(nspp_unknown),
            rel_cover_exotic = mean(rel_cover_exotic)) %>%
  ungroup() %>%
  dplyr::select(site, nspp_exotic, rel_cover_exotic,nspp_unknown) %>%
  pivot_longer(cols = c(nspp_exotic, rel_cover_exotic,nspp_unknown)) %>%
  ggplot(aes(x= value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free") +
  theme_classic()

all %>%
  ggplot(aes(x = nspp_unknown)) +
  geom_histogram() +
  theme_classic()
