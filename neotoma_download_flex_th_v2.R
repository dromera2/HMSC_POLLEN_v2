library(tidyverse)
library(terra)
library(neotoma2)
library(sf)

mask <- rast("./input/1140_bio.tif") %>% 
  subset(1)

mask_sf <- mask %>% 
  ext() %>%
  as.polygons() %>%
  st_as_sf()


st_crs(mask_sf) <- crs(mask)

records <- get_datasets(loc = mask_sf, datasettype = "pollen", all_data = TRUE) %>%
  neotoma2::filter(age_range_young >= 1) %>%
  get_downloads(all_data = TRUE)

count_by_site_NISP <- samples(records) %>%
  dplyr::filter(units == 'NISP' & 
                  elementtype == 'pollen' &
                  (ecologicalgroup %in% c('TRSH', 'UPHE', 'SEED', 'PALM', 'PLNT', 'UNID')) &
                  age <= 22000) 

myspecies <- c("Abies|A. pinsapo|Abies alba|A. alba|Cedrus atlantica|C. atlantica|Cedrus|Pinus halepensis|P. halepensis|Pinus nigra|P. nigra|Pinus pinea|P. pinea|Pinus sylvestris|P. sylvestris|evergreen|deciduous|Quercus faginea|Q. faginea|Quercus pyrenaica|Q. pyrenaica|Quercus ilex|Q. ilex|Quercus rotundifolia|Q. rotundifolia|Quercus coccifera|Q. coccifera")

taxa_names <- count_by_site_NISP %>%
  dplyr::filter(stringr::str_detect(variablename, myspecies)) %>%
  pull(variablename) %>% unique()

count_by_site_NISP <- count_by_site_NISP %>%
  mutate(variablename = replace(variablename, 
                                stringr::str_detect(variablename, "Abies|A. pinsapo|Abies alba|A. alba"),
                                "Abies"),
         variablename = replace(variablename,
                                stringr::str_detect(variablename, "Cedrus atlantica|C. atlantica|Cedrus"),
                                "Cedrus"),
         variablename = replace(variablename,
                                stringr::str_detect(variablename, "Pinus halepensis|P. halepensis|Pinus pinea|P. pinea"),
                                "Pinus halepensis/pinea type"),
         variablename = replace(variablename,
                                stringr::str_detect(variablename, "Pinus nigra|P. nigra|Pinus sylvestris|P. sylvestris"),
                                "Pinus nigra/sylvestris type"),
         variablename = replace(variablename,
                                stringr::str_detect(variablename, "Quercus ilex|Q. ilex|Quercus rotundifolia|Q. rotundifolia|Quercus coccifera|Q. coccifera|evergreen"),
                                "Quercus ilex/rotundifolia/coccifera type"),
         variablename = replace(variablename,
                                stringr::str_detect(variablename, "Quercus faginea|Q. faginea|Quercus pyrenaica|Q. pyrenaica|deciduous"),
                                "Quercus faginea/pyrenaica type")) %>%
  group_by(datasetid, age, variablename, units, ecologicalgroup, sitename, lat, long) %>%
  summarise(value = sum(value)) %>% as.data.frame()



#####remove sitenames without any specie of interest. currently disabled
sp <- c("Abies", "Cedrus",  "Pinus halepensis/pinea type", "Pinus nigra/sylvestris type", 
        "Quercus ilex/rotundifolia/coccifera type", "Quercus faginea/pyrenaica type")

# sites2keep <- count_by_site %>% dplyr::filter(variablename %in% sp) %>% pull(sitename) %>% unique()
# 
# count_by_site <- count_by_site %>% dplyr::filter(sitename %in% sites2keep)
####

breaks <- seq(1, 22000, 100)

count_by_site_NISP <- count_by_site_NISP %>%
  mutate(timebins = cut(count_by_site_NISP$age,
                        c(breaks, Inf),
                        labels = breaks, 
                        include.lowest = TRUE)) %>% 
  group_by(datasetid, sitename, long, lat, variablename, timebins) %>% 
  summarize(value = sum(value)) %>%
  ungroup()  %>% 
  arrange(datasetid, timebins)


totalcount_by_site_and_timebin <- count_by_site_NISP %>%
  group_by(datasetid, timebins) %>%
  summarize(totalpollen = sum(value)) %>%
  ungroup()

count_by_site_NISP <- count_by_site_NISP %>%
  left_join(totalcount_by_site_and_timebin) %>%
  mutate(polprop = round(100 * value / totalpollen, 1))

count_by_site_NISP <- count_by_site_NISP %>%
  pivot_wider(names_from = variablename, 
              values_from = polprop, 
              values_fill = 0) %>% 
  select(timebins, long, lat, sp) %>%
  na.omit()

thresholds <- count_by_site_NISP[,sp] %>% apply(2, FUN = function(x){
  max(x) * 0.05
})

count_by_site_NISP[,sp] <- count_by_site_NISP[,sp] %>% apply(2, FUN = function(x){
  total <- max(x)
  x = round(x * 100 / total, digits = 1)
  
  ifelse(x >= 5, 1, 0)
})


write.csv(count_by_site_NISP, "./input/neotoma_v2.csv")



