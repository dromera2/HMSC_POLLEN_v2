library(tidyverse)
library(terra)
library(usdm)
library(corrplot)

bioclimate <- rast("./input/1140_bio.tif") 
envirem <- rast("./input/1140_envirem.tif")

pH <- list.files(path = "./soil/1km/phh2o/", pattern = ".tif", full.names = TRUE )
sand <- list.files(path = "./soil/1km/sand/", pattern = ".tif", full.names = TRUE )
mask <- rast("./input/1140_bio.tif") %>% subset(1)


get.soil.data <- function(pH_list, sand_list, template) {
  
  pH <- rast(pH_list) %>% subset(c("phh2o_60-100cm_mean_1000",
                                   "phh2o_100-200cm_mean_1000"),
                                 negate = TRUE)
  
  sand <- rast(sand_list) %>% subset(c("sand_60-100cm_mean_1000",
                                       "sand_100-200cm_mean_1000"),
                                     negate = TRUE)
  
  # template <- template %>% 
  #   terra::unwrap() %>% 
  #   disagg(fact = 5, method = "bilinear")
  
  pH <- pH %>% 
    aggregate(fact = 10, fun = "mean", na.rm = TRUE) %>% 
    project(template, method = "bilinear")  %>% 
    app(mean)
  
  names(pH) <- "pH"
  
  sand <- sand %>% 
    aggregate(fact = 10, fun = "mean", na.rm = TRUE) %>% 
    project(template, method = "bilinear")  %>% 
    app(mean)
  
  names(sand) <- "sand"
  
  soil_data <- c(pH, sand)
  
  soil_data %>% scale()
  
  soil_data
}

soil <- get.soil.data(pH, sand, mask)


variables <- c(bioclimate, envirem, soil)

subset_vars <- variables %>%
  subset(c("bio2", "bio3", "bio4", "bio6", "bio7",
           "bio8", "bio9", "bio10", "bio11", "bio12", "bio13",
           "bio15", "bio16", "bio18", "bio19"),
         negate = TRUE)

names(subset_vars)

set.seed(1)

vif(subset_vars)

v1 <- vifcor(subset_vars, th = 0.7)
v2 <- vifstep(subset_vars, th = 10, method = "spearman")

v1 
v2 

env_stack <- subset_vars |> 
  subset(v2@excluded, negate = TRUE) 
plot(env_stack)

stats <- layerCor(env_stack, "cor", use = "complete.obs")
corr_matrix <- stats$correlation


jpeg("./figures/vars_correlation.jpg", width = 2000, height = 2000, res = 300)
corr_plot <- corrplot(corr_matrix, method = "ellipse")
dev.off()
