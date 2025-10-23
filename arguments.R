
myspecies <- c("Abies pinsapo", "Abies alba", "Cedrus atlantica", "Pinus halepensis", "Pinus pinea",
               "Pinus nigra", "Pinus sylvestris", "Quercus faginea", "Quercus pyrenaica", 
               "Quercus ilex", "Quercus rotundifolia", "Quercus coccifera")



lon_min <- -11.0416985
lon_max <- 12.0415425
lat_min <- 27.9583215
lat_max <- 44.0415905 


conf <- list(Latitude = "decimalLatitude",
             Longitude = "decimalLongitude",
             Date_collected = "eventDate",
             Scientific_name = "species")

mask <- rast("./input/1140_bio.tif") %>% subset(1)
bioclimate <- rast("./input/1140_bio.tif") 
envirem <- rast("./input/1140_envirem.tif")

var_list <- c("bio14", "aridityIndex", "continentality", "petSeasonality", "thermicityIndex")

#it is not needed to run this code, the data code will be provided.
# install.packages("usethis")
# usethis::edit_r_environ()
# 
# GBIF_USER="example_user"
# GBIF_PWD="example_password"
# GBIF_EMAIL="example_mail"
# 
# get.code <- function(myspecies) {
#   gbif_taxon_keys <- myspecies %>%
#     name_backbone_checklist() %>% # match to backbone
#     dplyr::pull(usageKey)
# 
# 
#   occurrence_data <- occ_download(
#     pred_in("taxonKey", gbif_taxon_keys),
#     #### pred("basisOfRecord", "PRESERVED SPECIMEN"),
#     pred("hasCoordinate", TRUE),
#     pred("hasGeospatialIssue", FALSE),
#     format = "SIMPLE_CSV"
#   )
# 
#   occurrence_data
# }
# 
# get.code(myspecies)
# 
# Download this data file from:  https://www.gbif.org/occurrence/download/0113662-240626123714530
# code <- '0113662-240626123714530'

# my_download_metadata <- occ_download_meta(code)
# gbif_citation(my_download_metadata)

# gbif <- get.gbif.data(code, minlon = lon_min, maxlon = lon_max,
#                          minlat = lat_min, maxlat = lat_max)
# 
# write.csv(gbif, "./input/gbif.csv")
gbif_names <- read.csv("./input/gbif.csv") %>% select(-1) %>% colnames()

gbif_data <- read_sf("./input/gbif_data_v2.shp") %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  select(-1)

colnames(gbif_data) <- gbif_names


prev <- gbif_data %>%
  select(decimalLongitude, decimalLatitude, species) %>% 
  mutate(value = 1) %>% 
  pivot_wider(names_from = species, values_from = value, values_fill = 0) %>% 
  select(myspecies) %>% 
  apply(2, mean) %>% 
  as.data.frame() %>% 
  mutate(sp = myspecies) 

colnames(prev) <- c("spec_sens", "sp")#Note that column 1 is the prevalence, we just called it spec_sens to ease the reuse of functions

prev_th <- list()
prev_th[[1]] <- prev
prev_th[[2]] <- prev


# Download data from this source: https://files.isric.org/soilgrids/former/2017-03-10/aggregated/1km/
pH <- list.files(path = "./soil/1km/phh2o/", pattern = ".tif", full.names = TRUE )
sand <- list.files(path = "./soil/1km/sand/", pattern = ".tif", full.names = TRUE )
  
i <- 1:10
# i <- 1:2


samples <- 100
nchains <- 4
thin <- 50
transient <- 80*thin
verbose = 100
nParallel = 2


coord_names <- c("x", "y")
var_names <- c(var_list, "pH", "sand")

formula <- ~poly(bio14, degree = 2, raw = TRUE) + 
  poly(aridityIndex, degree = 2, raw = TRUE) +
  poly(continentality, degree = 2, raw = TRUE) +
  poly(petSeasonality, degree = 2, raw = TRUE) +
  poly(thermicityIndex, degree = 2, raw = TRUE) +
  poly(pH, degree = 2, raw = TRUE) +
  poly(sand, degree = 2, raw = TRUE) 

groupnames = c("climate", "soil")

group <- c(rep(1, 11), rep(2, 4))


pollen_sites <- read.csv("./input/moroccotrees.csv") |>
  dplyr::select(1, 3, 4)

colnames(pollen_sites) <- c("species", "decimalLongitude", "decimalLatitude")

trees_data <- read.csv("./input/moroccotrees.csv")
pollen_data <- read.csv("./input/MoroccoPollenData.csv")

starts <- seq(-22000, -1000, 1000)
ends <- seq(-21001, -1, 1000)


y_values <- seq(-22000, -1, 100)

rcps <- c("rcp2.6", "rcp4.5", "rcp6.0", "rcp8.5")
gcms <- c("CESM1-CAM5", "CSIRO-Mk3-6-0", "IPSL-CM5A-MR")

colors <- c("#009E73", rep("#F0E442", 2))

records <- read.xlsx("./input/REGISTROS.xlsx", sheetName = "complete") |> 
  select(-site)

neotoma <- read.csv("./input/neotoma_v2.csv") %>%
  select(-1) |> 
  bind_rows(records)

neotoma_2_eval <- read.csv("./input/neotoma_v2.csv", check.names = FALSE) %>% select(-1)

colnames(records) <- colnames(neotoma_2_eval)

neotoma_2_eval <- neotoma_2_eval |> 
  bind_rows(records)
