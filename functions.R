

# remove records without coordinates
clean.coords <- function(myspecies_coords) {
  myspecies_coords <- myspecies_coords %>%
    filter(!is.na(decimalLongitude)) %>%
    filter(!is.na(decimalLatitude)) %>% 
    mutate(countryCode = replace(countryCode, countryCode == "", "SK"))
  
  
  
  #convert country code from ISO2c to ISO3c
  myspecies_coords$countryCode <-  countrycode(myspecies_coords$countryCode, 
                                               origin =  'iso2c',
                                               destination = 'iso3c',
                                               custom_match = c("XK" = "SRB"))
  
  flags <- clean_coordinates(x = myspecies_coords, 
                             lon = "decimalLongitude", 
                             lat = "decimalLatitude",
                             # countries = "countryCode",
                             species = "species",
                             tests = c("capitals", "centroids", "equal","gbif",
                                       "institutions",
                                       # "outliers", 
                                       # "seas", 
                                       "zeros", 
                                       # "urban", 
                                       "validity", 
                                       "duplicates")) # most test are on by default
  
  #Exclude problematic records
  dat_cl <- myspecies_coords[flags$.summary,]
  dat_cl
}



#Remove records with low coordinate precision
low.precision <- function(dat_cl) {
  dat_cl <- dat_cl %>%
    filter(coordinateUncertaintyInMeters / 1000 <= 10 | is.na(coordinateUncertaintyInMeters))
  dat_cl
}


individual.count <- function(dat_cl) {
  dat_cl <- dat_cl%>%
    filter(individualCount > 0 | is.na(individualCount))%>%
    filter(individualCount < 99 | is.na(individualCount))
  dat_cl
}


record.age <- function(dat_cl) {
  dat_cl <- dat_cl%>%
    filter(year > 1945) # remove records from before second world war
  dat_cl
}

long_lat.cut <- function(dat_cl, minlon, maxlon, minlat, maxlat) {
  dat_cl <- dat_cl%>%
    filter (decimalLongitude > minlon) %>% 
    filter(decimalLongitude < maxlon) %>% 
    filter(decimalLatitude > minlat) %>% 
    filter(decimalLatitude < maxlat)
  dat_cl
}


get.gbif.data <- function(code, minlon, maxlon, minlat, maxlat) {
  
  occurrences <- occ_download_get(code) %>%
    occ_download_import()
  
  filtered_longlat <- long_lat.cut(occurrences, minlon, maxlon, minlat, maxlat)
  filtered_precision <- low.precision(filtered_longlat)
  filtered_counts <- individual.count(filtered_precision)
  filtered_coordinates <- clean.coords(filtered_counts)

  filtered_coordinates
  
}

make.mask <- function(indf, resolution) {
  r <- terra::rast()
  ext(r) <- c(min(indf$Longitude) %>% floor(),
                   max(indf$Longitude) %>% ceiling(),
                   min(indf$Latitude) %>% floor(),
                   max(indf$Latitude) %>% ceiling())
  res(r) <- resolution
  crs(r) <- "epsg:4326"
  
  terra::wrap(r)
}
########blockcv

get.blocks <- function(data, n = 3, i = 100, size = 50000) {
  data <- data %>% select(x, y) 
  data$occ <- 1
  
  sf_data <- st_as_sf(x = data, 
                      coords = c("x", "y"),
                      crs = "epsg:4326")
  
  #Spatial blocking by specified range with random assignment
  set.seed(1)
  sb <- cv_spatial(x = sf_data,
                   column = "occ",
                   size = size,
                   k = n,
                   selection = "random",
                   iteration = i, # find evenly dispersed folds
  )
  blocks <- data.frame(partition = sb$folds_ids)
  
  blocks
}


##############

completeness <- function (indf, raster, recs = 10) {
  require(tidyverse)
  require(terra)  
  
  raster <- terra::unwrap(raster)
  indf <- as.data.frame(indf)
  
  coords <- indf %>% 
    select(Longitude, Latitude) 
  
  indf$cust_cell_id <- cellFromXY(raster, coords)
  
  dat1 <- indf %>% # Remove duplicated entries
    select(Scientific_name, Date_collected, cust_cell_id) %>% 
    distinct(Scientific_name, Date_collected, cust_cell_id)
  dat2 <- dat1 %>% group_by(cust_cell_id) %>% # Number of records in each gridcell
    summarize(cell_ct = n())
  dat3 <- dat2 %>% # Filter by minimum  number of records
    filter(cell_ct >= recs) 
  if (nrow(dat3) < 1) {
    stop("Too few data records to compute completeness")
  }
  retmat <- NULL
  for (i in 1:nrow(dat3)) {
    print(paste0(i, " / ", nrow(dat3)))
    cust_cell_id <- dat3$cust_cell_id[i] # The cell ID
    nrec <- dat3$cell_ct[i] # Number of records in the cell
    cset <- dat1[which(dat1$cust_cell_id == cust_cell_id), ] # Comm. matrix of the cell
    csum <- cset %>% group_by(Scientific_name) %>% summarize(sct = n()) # Number of records for each species in the cell
    Q1 <- csum %>% filter(sct == 1) %>% nrow() # Number of unique species
    Q2 <- csum %>% filter(sct == 2) %>% nrow() # Number of duplicate species
    m <- cset %>% distinct(Date_collected) %>% nrow() # Number of samples (times sampled the grid) 
    Sobs <- cset %>% distinct(Scientific_name) %>% nrow()
    if (Sobs > 0) {
      Sest <- Sobs + (((m - 1)/m) * ((Q1 * (Q1 - 1))/(2 * (Q2 + 1))))
      c <- Sobs/Sest
      retmat <- rbind(retmat, c(cust_cell_id, nrec, Q1, Q2, m, Sobs, Sest, c))
    }
  }
  retmat <- as.data.frame(retmat)
  names(retmat) <- c("Cell_id", "nrec", "Q1", "Q2", "m", "Sobs", "Sest", "c")
  
  retmat
}

rasterize_completeness <- function(comp_an, var, raster) {
  
  raster <- raster %>% terra::unwrap()
  
  new_values <- rep(NA, ncell(raster))
  new_values[comp_an$Cell_id] <- comp_an[, var]
  values(raster) <- new_values
  raster %>% terra::wrap()
}

sample.data <- function(complection_rast, gbif_data, complection_data,  threshold) {
  raster <- complection_rast %>% terra::unwrap()
  
  filtered_data <- gbif_data %>% as.data.frame() 
  
  coords <- filtered_data %>% select(decimalLongitude, decimalLatitude)
  
  filtered_data$Cell_id <- cellFromXY(raster, coords)
  
    comp_data <- complection_data %>%
      filter(c >= threshold) %>% 
      select(Cell_id) 
  
 
    filtered_th <- filtered_data %>%
      filter(!(species %in% c("Abies pinsapo", "Cedrus atlantica"))) %>% 
      filter(Cell_id %in% comp_data$Cell_id)
    
    filtered_under <- filtered_data %>% 
      filter(species %in% c("Abies pinsapo", "Cedrus atlantica"))
    
  filtered_data <- bind_rows(filtered_th, filtered_under)
  
  weight_sample <- complection_rast %>% weight.rast(threshold)
  
  output <- left_join(filtered_data, weight_sample, by = "Cell_id") %>% 
    mutate(weight = replace(weight, is.na(weight), 0))
  
  output
  
}

under.threshold <- function(raw_data, filtered_data){
  under_data <- raw_data %>% 
    filter(!(gbifID %in% filtered_data$gbifID))
  
  under_data
}

get.climate.data <- function(bio, env, variables) {
  
  clim <- c(bio, env) |> subset(variables) |> terra::wrap()
  clim
}


get.soil.data <- function(pH_list, sand_list, template) {
  
  pH <- rast(pH_list) %>% subset(c("phh2o_60-100cm_mean_1000",
                                   "phh2o_100-200cm_mean_1000"),
                                 negate = TRUE)
  
  sand <- rast(sand_list) %>% subset(c("sand_60-100cm_mean_1000",
                                       "sand_100-200cm_mean_1000"),
                                     negate = TRUE)
  
  
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
  
  soil_data %>% scale() %>%  terra::wrap()
}

na2idw <- function(varname, data, template, formula) {
  mg <- gstat(id = varname, formula = formula, locations = ~x+y, data=as.data.frame(data, xy = TRUE), 
              nmax=7, set=list(idp = .5))
  z <- terra::interpolate(data, mg, debug.level=0, index=1, xyNames = c("x", "y"))
  z<- terra::cover(data, z)
  z <- mask(z, template)
  z
}

get.Suse.data <- function(var, template, name, formula) {
  var <- var %>% aggregate(fact = 10) %>% project(template, method = "bilinear")
  names(var) <- name
  
  var <- na2idw(name, var, template, formula)
  var %>% scale() %>% terra::wrap()
}


filter.data <- function(type, data, comp_rast, comp_data, th){
  if(type %in% c("over", "under")){
    source_data <- data
    data <- sample.data(comp_rast, data, comp_data, th)
    if (type == "under"){
      data <- under.threshold(source_data, data)
    }
  }
  return(data)
}



prepare.data <- function(type, occ_data, comp_rast, comp_data, th, clim_data, soil_data,
                         cut, sample = NULL, i, species_list) {
  
  
  base_data <- occ_data
  
  occ_data <- filter.data(type, occ_data, comp_rast, comp_data, th)
  
  points <- occ_data %>% dplyr::select(c(species, decimalLongitude, decimalLatitude, weight))
  
  points <- points  %>%  
    mutate(weight = coalesce(weight, 0)) %>% 
    mutate(species = na_if(species, "")) %>% 
    na.omit() %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"))
  
  clim_data <- clim_data %>% terra::unwrap()
  
  template <- clim_data %>% subset(1)
  
  points <- points %>% mutate(pixelId = cellFromXY(template, points %>% st_coordinates()), occurrence = 1) %>% 
    dplyr::select(species, pixelId, occurrence, weight) %>% st_drop_geometry() %>% 
    pivot_wider(names_from = "species", values_from = "occurrence", values_fn = max, values_fill = 0) %>% 
    distinct(pixelId, .keep_all = TRUE)
  
  coords <- template %>% xyFromCell(cell = points$pixelId) %>% as.data.frame()
  
  clim <- clim_data %>%  terra::extract(y = points$pixelId)
  
  soil <- soil_data %>% unwrap() %>% terra::extract(y = points$pixelId)
  
  model_data <- bind_cols(coords, clim, soil, points |> relocate(c("pixelId", "weight", species_list))) %>% 
    dplyr::select(-pixelId) %>% na.omit()
  
  if(cut == TRUE){
    set.seed(i)  
    
    low_prevalent <- model_data %>% filter(`Abies pinsapo` == 1 | `Cedrus atlantica` == 1)
    
    subset_data <- model_data %>%
      filter(`Abies pinsapo` == 0) %>% 
      filter(`Cedrus atlantica` == 0)
    
    part_sample <- sample - nrow(low_prevalent)
    
    subset_data <- subset_data %>% slice_sample(n = part_sample, weight_by = subset_data$weight) 
    
    model_data <- bind_rows(low_prevalent, subset_data)
    
    
  }
  
  points <- base_data %>% dplyr::select(c(species, decimalLongitude, decimalLatitude))
  
  points <- points  %>%  
    mutate(species = na_if(species, "")) %>% 
    na.omit() %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"))
  
  points <- points %>% mutate(pixelId = cellFromXY(template, points %>% st_coordinates()), occurrence = 1) %>% 
    dplyr::select(species, pixelId, occurrence) %>% st_drop_geometry() %>% 
    pivot_wider(names_from = "species", values_from = "occurrence", values_fn = max, values_fill = 0) %>% 
    distinct(pixelId, .keep_all = TRUE)
  
  coords <- template %>% xyFromCell(cell = points$pixelId) %>% as.data.frame()
  
  clim <- clim_data %>%  terra::extract(y = points$pixelId)
  
  soil <- soil_data %>% unwrap() %>% terra::extract(y = points$pixelId)
  
  base_data <- bind_cols(coords, clim, soil, points |> relocate(c("pixelId", species_list))) %>% 
    dplyr::select(-pixelId) %>% na.omit() %>% add_column(weight = 0, .after = "sand")
  
  anti_data <- anti_join(base_data, model_data, by = c("x", "y"))
  
  random_points <- anti_data %>% slice_sample(n = nrow(model_data)/4)
  
  model_data <- bind_rows(model_data, random_points) %>% na.omit()
  
  mask_rast <- ifel(!is.na(template), 1, 0)
  
  presences <- base_data %>% select(x, y) %>% st_as_sf(coords = c("x" ,"y"), crs = crs(template))
  reproj_presences <- st_transform(presences, crs = "+proj=aeqd +lat_0=41.5 +lon_0=-73.5 +datum=WGS84")
  
  buffer <- st_buffer(reproj_presences, dist = 50000) %>%
    st_transform(crs = crs(template))
  
  mask_pol <- mask_rast %>%
    as.polygons() %>%
    st_as_sf() %>%
    st_transform(crs = crs(template)) %>%
    select(geometry) %>%
    slice(2)
  
  set.seed(i)
  absences_points <- st_sample(st_as_sfc(st_bbox(mask_rast)), size = nrow(model_data), type = "random") %>%
    st_as_sf() %>%
    st_intersection(mask_pol) %>%
    st_difference(st_union(buffer)) %>% 
    slice_sample(n = nrow(model_data)/5)
  
  absences_coords <- absences_points %>% st_coordinates() %>% as.data.frame()
  
  absences_cells <- absences_points %>% mutate(pixelId = cellFromXY(template, absences_points %>% st_coordinates()))
  
  absences_clim <- clim_data %>%  terra::extract(y = absences_cells$pixelId)
  
  absences_soil <- soil_data %>% unwrap() %>% terra::extract(y = absences_cells$pixelId)
  
  absences_data <- bind_cols(absences_coords, absences_clim, absences_soil)
  
  data_shell <- model_data %>%
    slice_sample(n = nrow(absences_data))
  
  data_shell[,1:9] <- absences_data
  data_shell[,10:22] <- 0
  
  model_data <- bind_rows(model_data, data_shell)
  model_data <- model_data %>%
    distinct(x, y, .keep_all = TRUE) %>%
    na.omit()
  
  model_data
  
  partition <- get.blocks(model_data)
  
  model_data$partition <- partition$partition
  
  model_data
}


prepare.env.data <- function(start, end){
  
  
  f <- function(x, format = "%B") { 
    months = format(as.Date(paste0("01-", 1:12, "-1970")), format)
    factor(format(x, format), levels = months)
  }
  
  time.aggregate <- function(var, time) {
    
    temp <- list()
    
    for(i in time) {
      temp[[match(i, time)]] <- var |> 
        slice("time", i:(i+1199)) |> 
        aggregate(f, mean) |> 
        st_set_dimensions(1, values = lubridate::ymd(paste0(sprintf("%04d", ( (i-1) / 120 ) + 1), "-", 1:12, "-", 1)), names = "time")
    }
    temp
  }
  
  
  time <- seq (1, (end - start + 1)*12 , 1200)
  
  
  tmax <- dsclimtools::read_dsclim("./dsclim", "tasmax", start, end, calendar_dates = TRUE, proxy = FALSE)
  
  tmax_agg <- tmax |> 
    time.aggregate(time) |>
    lapply(FUN = function(x){x |> rast() |> stack()}) |> 
    lapply(FUN = function(x){
      crs(x) <- CRS('+init=EPSG:4326')
      x
    })
  
  rm(tmax)
  gc()
  
  tmin <- dsclimtools::read_dsclim("./dsclim", "tasmin", start, end, calendar_dates = TRUE, proxy = FALSE) 
  
  tmin_agg <- tmin |>
    time.aggregate(time) |> 
    lapply(FUN = function(x){x |> rast() |> stack()}) |> 
    lapply(FUN = function(x){
      crs(x) <- CRS('+init=EPSG:4326')
      x
    })
  
  rm(tmin)
  gc()
  
  t <- dsclimtools::read_dsclim("./dsclim", "tas", start, end, calendar_dates = TRUE, proxy = FALSE)
  
  t_agg <- t |>
    time.aggregate(time) |>
    lapply(FUN = function(x){x |> rast() |> stack()}) |>
    lapply(FUN = function(x){
      crs(x) <- CRS('+init=EPSG:4326')
      x
    })
  
  rm(t)
  gc()
  
  prec <- dsclimtools::read_dsclim("./dsclim", "pr", start, end, calendar_dates = TRUE, proxy = FALSE)
  
  prec_agg <- prec |>
    time.aggregate(time) |>
    lapply(FUN = function(x){x |> rast() |> stack()}) |> 
    lapply(FUN = function(x){
      crs(x) <- CRS('+init=EPSG:4326')
      x
    })
  
  rm(prec)
  gc()
  
  bio_list<- list()
  
  for(t in 1:length(prec_agg)){
    
    bio_list[[t]] <- dismo::biovars(prec_agg[[t]], tmin_agg[[t]], tmax_agg[[t]]) |> rast()
  }
  etsolrad_seq <- seq(start, end, 10)
  
  etsolrad <- lapply(etsolrad_seq,
                     FUN = function(x, y){envirem::ETsolradRasters(y, year = x)},
                     prec_agg[[1]] |> subset(1) |> rast())
  
  etsolrad <- do.call(mean, etsolrad)
  
  clim_list <- list()
  
  for(c in 1:length(prec_agg)){
    
    prec <- prec_agg[[c]] |> rast()
    terra::set.names(prec, paste0("precip_", 1:12))
    
    t <- t_agg[[c]] |> rast()
    terra::set.names(t, paste0("tmean_", 1:12))
    
    tmin <- tmin_agg[[c]] |> rast()
    terra::set.names(tmin, paste0("tmin_", 1:12))
    
    tmax <- tmax_agg[[c]] |> rast()
    terra::set.names(tmax, paste0("tmax_", 1:12))
    
    clim_list[[c]] <- c(prec, t, tmin, tmax)
  }
  
  envirem_data <- lapply(clim_list, FUN = function(x){envirem::generateEnvirem(x, etsolrad, var = "all")}) |> 
    lapply(FUN = function(x){terra::wrap(x)})
  
  bio_list <- bio_list |> lapply(FUN = function(x){terra::wrap(x)})
  
  
  variables <- list()
  variables[[1]] <- bio_list
  variables[[2]] <- envirem_data
  variables
}

#######


prepare.site.data <- function(points, clim_data, soil_data) {
  
  points <- points  |> 
    st_as_sf(coords = c("long", "lat"))
  
  clim_data <- clim_data %>% terra::unwrap()
  
  template <- clim_data %>% subset(1)
  
  points <- points %>% mutate(pixelId = cellFromXY(template, points %>% st_coordinates())) %>% 
    dplyr::select(pixelId) |>  
    distinct(pixelId, .keep_all = TRUE)
  
  coords <- template %>% xyFromCell(cell = points$pixelId) %>% as.data.frame()
  
  clim_data <- clim_data %>%  terra::extract(y = points$pixelId)
  
  soil_data <- soil_data %>% unwrap() %>% terra::extract(y = points$pixelId)
  
  model_data <- bind_cols(coords, clim_data, soil_data) |>  na.omit()
  
  model_data
}



prepare.timeseries.data <- function(points, clim_data, soil_data) {
  
  clim_data <- do.call(Map, c(c, clim_data))
  
  points <- points  |> 
    st_as_sf(coords = c("long", "lat"))
  
  template <- clim_data[[1]][[1]] |> terra::unwrap() |>  subset(1)
 
  timeseries_list <- list()
  
  for(i in 1:length(clim_data[[1]])){
    
    bio <- clim_data[[1]][[i]] |> terra::unwrap() |> subset("bio14")
    env <- clim_data[[2]][[i]] |> terra::unwrap() |>  subset(c("aridityIndexThornthwaite", "continentality", "PETseasonality", "thermicityIndex"))
    names(env) <- c("aridityIndex", "continentality", "petSeasonality", "thermicityIndex")
    
    timeseries_list[[i]] <- c(bio, env)
  }
  
  
  points <- points %>% mutate(pixelId = cellFromXY(template, points %>% st_coordinates())) %>% 
    dplyr::select(pixelId) |>  
    distinct(pixelId, .keep_all = TRUE)
  
  coords <- template %>% xyFromCell(cell = points$pixelId) %>% as.data.frame()
  
 
  clim <- timeseries_list |> lapply(function(x){terra::extract(x, y = points$pixelId)})
  
  soil_data <- soil_data %>% unwrap() %>% terra::extract(y = points$pixelId)
 
  soil_list <- list()
  coords_list <- list()
  
  for (s in 1:length(clim)){
    soil_list[[s]] <- soil_data
    coords_list[[s]] <- coords
  }
  
  model_data <- list()
  
  model_data[[1]] <- clim
  model_data[[2]] <- soil_list
  model_data[[3]] <- coords_list
  
  model_data
}



plot.diagnostics <- function(model, name, i) {
  mpost = convertToCodaObject(model)
  
  ess.beta = effectiveSize(mpost$Beta)
  
  ess.V = effectiveSize(mpost$V)
  
  psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  
  psrf.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
  
  par(mfrow = c(2, 2))
  pdf(paste0("figures/diagnostics/diagnostics_", name, i, ".pdf"))
  hist(ess.beta, xlab = expression("Effective sample size" ~ beta ~ ""), main=NULL, col= "light green")
  hist(ess.V, xlab = expression("Effective sample size" ~ v ~ ""))
  hist(psrf.beta, xlab = expression("Potential scale reduction factor" ~ beta ~ ""), main=NULL, col= "light green")
  hist(psrf.V, xlab = expression("Potential scale reduction factor" ~ v ~ ""))
  dev.off()
  
}


plot.variance <- function(var, i, cols, name) {
  
  temp <- list()
  
  for(d in 1:length(i)){
    values <- var[[d]]$vals
    
    temp[[d]] <- values
  }
  
  mean <- (Reduce(`+`, temp) / length(temp)) %>% 
    t() %>%  as.data.frame() %>% 
    mutate(species = rownames(.)) %>% 
    pivot_longer(cols = cols) 
  
  
  if(name == "grouped"){
    
    mean <- mean %>% 
      mutate(name = case_when(
        name == "climate" ~ "Climate",
        name == "soil" ~ "Soil",
        .default = name))
  } else{
    mean
  }
  
  mean %>% ggplot() +
    geom_col(aes(x = value, y = species, fill = name)) +
    labs(x = "Explained variance (proportion)", y = "") +
    theme(
      axis.text = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      legend.title = element_blank(),
      legend.text = element_text(size = 15)
    ) 
  
  if(name == "grouped"){
    ggsave(paste0("./figures/variance_partitioning_", name, "_V3.jpg"), height = 4, width = 8)
  } else{
    ggsave(paste0("./figures/variance_partitioning_", name, "_V3.jpg"), height = 4, width = 12)
    
  }
  
  mean
  
}


spatial.prediction <- function(climate, soil, species, model, i, dataset) {

  climate <- climate %>% terra::unwrap() #%>% crop(mask)
  soil <- soil %>% terra::unwrap() #%>% crop(mask)
  
  env_data <- c(climate, soil)
  
  XData.grid <- env_data %>% terra::values(dataframe = TRUE)
  xy.grid <- env_data %>% xyFromCell(cell = 1:ncell(.)) %>% as.data.frame()
  
  complete.index <- which(complete.cases(XData.grid))
  Gradient <- prepareGradient(model, XDataNew = XData.grid[complete.index,], sDataNew = list(units=xy.grid[complete.index,]))
  
  predY <- predict(model, Gradient = Gradient, expected = TRUE, predictEtaMean = TRUE, nParallel=2)
  
  EpredY <- matrix(NA, nrow(XData.grid), length(species))
  
  library(abind)
  EpredY[complete.index,] <- apply(abind(predY,along=3),c(1,2),mean)
  
  mapData <- data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY)
  colnames(mapData) <- c("x", "y", species)
  library(reshape2)
  mapData <- melt(mapData, id.vars = c("x", "y"))
  
  mapData$i <- i
  mapData$dataset <- dataset
  
  mapData
 
}

coords.prediction <- function(comm, vars, coords, species, model, i, dataset) {
  
  XData.grid <- comm %>% select(vars)
  xy.grid <- comm %>% select(coords)
  
  complete.index <- which(complete.cases(XData.grid))
  Gradient <- prepareGradient(model, XDataNew = XData.grid[complete.index,], sDataNew = list(units=xy.grid[complete.index,]))
  
  predY <- predict(model, Gradient = Gradient, expected = TRUE, predictEtaMean = TRUE, nParallel=2)
  
  EpredY <- matrix(NA, nrow(XData.grid), length(species))
  
  library(abind)
  EpredY[complete.index,] <- apply(abind(predY,along=3),c(1,2),mean)
  
  mapData <- data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY)
  colnames(mapData) <- c("x", "y", species)
  library(reshape2)
  mapData <- melt(mapData, id.vars = c("x", "y"))
  
  mapData$i <- i
  mapData$dataset <- dataset
  
  mapData
  
}

timeseries.prediction <- function(comm, vars, coords, species, model, i, dataset) {
  
  clim_list <- comm[[1]]
  soil_list <- comm[[2]]
  coords_list <- comm[[3]]
  
  data <- list()
  
  for(f in 1:length(clim_list)){
  data[[f]] <- bind_cols(clim_list[[f]], soil_list[[f]], coords_list[[f]])
  }
  
  prediction <- data |> lapply(FUN = function(x, y, z, a, b, c, w){
    coords.prediction(x, vars, coords, species, model, i, dataset)
  })
  
  prediction <- bind_rows(prediction, .id = "time")
  
}


prediction.map <- function(pred, sp_names, datasets, name) {
  
  
  pred <- pred %>% group_by(x, y, variable, dataset) %>% summarize(value = mean(value, na.rm = TRUE))
  
  
  pred %>% 
    ggplot(aes(x=x, y=y, fill= value, col = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_colour_viridis_c() +
    theme(legend.position="bottom", aspect.ratio= 1) +
    facet_wrap(. ~ variable) +
    theme_bw() +
    theme(legend.title=element_blank()) +
    coord_equal()
  
  ggsave(paste0("./figures/", name, ".jpg"), width = 10, height = 8)
  pred
}

get.betas <- function (model, species, type, i, vars){
  
  postBeta <- getPostEstimate(model, parName="Beta")
  beta.mod <- data.frame(t(postBeta$mean))
  colnames(beta.mod) <- c("int",  c(paste0(vars, "_", 1), paste0(vars, "_", 2))) 
  
  beta.mod <- melt(beta.mod, measure.vars= c("int",  c(paste0(vars, "_", 1), paste0(vars, "_", 2)))) #idem
  beta.mod$species <- species
  beta.mod$type <- type
  beta.mod$i <- i
  colnames(beta.mod) <- c("variable", "beta", "species", "type", "i")
  return(beta.mod)
}

get.omegas <- function(model, sp_order){
  
  OmegaCor = computeAssociations(model)
  
  plotOrder <- sp_order
  
  omega <- OmegaCor[[1]]$mean[plotOrder,plotOrder]
} 

omegas_to_df <- function(matrix, sp_order, type, i) {

  plotOrder <- sp_order
  
  prob <- matrix[plotOrder, plotOrder]
  prob_matrix <- as.matrix(prob)
  prob_matrix[upper.tri(prob_matrix, diag = TRUE)] <- NA
  prob_matrix <- na.omit(melt(prob_matrix))
  prob_matrix <- cbind(prob_matrix, paste0(prob_matrix[,1], "-", prob_matrix[,2]))
  prob_matrix[,5] <- type
  prob_matrix[,"i"] <- i
  colnames(prob_matrix) <- c("var1", "var2", "prob", "sp-sp", "type", "i")
  prob_matrix
}

plot.omegas <- function(corr_plot, sp_order) {
  
  corr_plot <- corr_plot %>% group_by(type, var1, var2) %>% summarize(prob = mean(prob))
  
  plotOrder <- sp_order
  
  corr_plot$var1 <- corr_plot$var1 %>% factor(levels = plotOrder)
  corr_plot$var2 <- corr_plot$var2 %>% factor(levels = rev(plotOrder))
  
  ggplot(corr_plot) +
    geom_raster(aes(x = var1, y = var2, fill = prob)) +
    scale_fill_gradientn(colors= c("#7f3b08", "#e08214", "#f5f5f5", "#8073ac", "#2d004b"),
                         values = scales::rescale(c(-1, -0.75, 0, 0.75, 1)), limits= c(-1, 1)) +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(legend.title= element_blank(), axis.text.x=element_text(angle = 45, hjust = 1)) +
    theme(text = element_text(size=rel(4.5)), legend.text = element_text(size=10),
          legend.position = "inside",
          legend.position.inside = c(0.1, 0.2))
  
  
  filename <- "./figures/omegas_V3.jpg"
  ggplot2::ggsave(filename, device = "jpg", width = 7, height = 6)
  
  corr_plot
}

plot.betas <- function(data){
  data <- data %>% 
    group_by(species, variable, type) %>% 
    summarize(beta = mean(beta)) %>% 
    filter(variable != "int")
  
  data %>% ggplot() +
    geom_tile(aes(x = factor(type, levels = c("under", "full", "over")), y = species, fill = beta)) +
    scale_fill_gradientn(colors= c("#8e0152", "#c51b7d", "#f5f5f5", "#4d9221", "#276419"), values = c(0, 0.4, 0.5, 0.6, 1), limits= c(-1, 1)) +
    facet_wrap(.~factor(variable, levels = c("PC1_1", "PC1_2", "PC2_1", "PC2_2", "PC3_1", "PC3_2",
                                             "sand_1", "sand_2", "pH_1", "pH_2"))) +
    xlab("") +
    theme_bw()
  
  ggsave("./figures/betas.pdf", width = 9, height = 6)
  
  data %>% 
    filter(variable %in% c("PC1_1", "PC1_2", "PC2_1", "PC2_2", "sand_1", "pH_1")) %>% 
    ggplot() +
    geom_tile(aes(x = factor(type, levels = c("under", "full", "over")), y = species, fill = beta)) +
    scale_fill_gradientn(colors= c("#8e0152", "#c51b7d", "#f5f5f5", "#4d9221", "#276419"), values = c(0, 0.4, 0.5, 0.6, 1), limits= c(-1, 1)) +
    facet_wrap(.~factor(variable, levels = c("PC1_1", "PC1_2", "PC2_1", "PC2_2", "sand_1", "pH_1"))) +
    xlab("") +
    theme_bw()
  
  ggsave("./figures/subset_betas.pdf", width = 9, height = 6)
  
  return(data)
}

##########



neotoma.validation <- function(neotoma, mask, timeseries_pred, CM = TRUE) {
  
  pollen_types <- neotoma %>% colnames() %>% .[seq_along(.) %in% c(4:9)]
  
  neo_points <- neotoma  %>%  
    pivot_longer(cols = pollen_types,
                 names_to = "variablename",
                 values_to = "value") %>% 
    st_as_sf(coords = c("long", "lat"))
  
  neo_points <- neo_points %>% mutate(pixelId = cellFromXY(mask, neo_points %>% st_coordinates())) 
  
  neotoma_data <- neo_points %>% 
    st_drop_geometry() %>% 
    group_by(pixelId, variablename, timebins) %>% 
    summarize(value = max(value)) %>% 
    ungroup() 
  
  if(CM != TRUE){
    timeseries_pred <- timeseries_pred %>% mutate(dataset = "over")
    
    colnames(timeseries_pred) <- c("time", "x", "y", "value", "i", "variable", "dataset")  
  }else{
    timeseries_pred
  }
  
  pred_points <- timeseries_pred   %>%  
    na.omit() %>%
    st_as_sf(coords = c("x", "y"))
  
  pred_points <- pred_points %>% mutate(pixelId = cellFromXY(mask, pred_points %>% st_coordinates())) %>% st_drop_geometry()
  
  pred_points <- pred_points %>% select(-dataset) %>% group_by(pixelId, variable, time) %>% 
    summarize(value = mean(value)) %>% ungroup() 
  
  ###no pred for this pixel
  
  pred_pixels <- pred_points %>% pull(pixelId) %>% unique()
  
  neotoma_pixels <- neotoma_data %>% pull(pixelId) %>% unique()
  
  pixel2keep <- intersect(pred_pixels, neotoma_pixels)
  
  neotoma_data <- neotoma_data %>% dplyr::filter(pixelId %in% pixel2keep)
  
  
  pred_points <- pred_points %>% dplyr::filter(pixelId %in% pixel2keep)
  
  
  pred_points <- pred_points %>% mutate(time = as.numeric(time))
  
  pred_points <- pred_points[order(pred_points$time),]
  
  sites <- pred_points %>% pull(pixelId) %>% unique()
  
  species <- pred_points %>% pull(variable) %>% unique() %>% as.character()
  
  
  val <- data.frame(age = seq(1, 22000, 100), time = seq(220, 1, -1))
  pred_points <- left_join(pred_points, val, by = "time") %>% select(-time)
  
 
  pred_points <- pred_points %>%
    mutate(variable = 
             case_when(variable == "Abies pinsapo" ~ "Abies",
                       variable == "Abies alba" ~ "Abies",
                       variable == "Cedrus atlantica" ~ "Cedrus",
                       variable == "Pinus halepensis" ~ "Pinus halepensis/pinea type",
                       variable == "Pinus pinea" ~ "Pinus halepensis/pinea type",
                       variable == "Pinus nigra" ~ "Pinus nigra/sylvestris type",
                       variable == "Pinus sylvestris" ~ "Pinus nigra/sylvestris type",
                       variable == "Quercus faginea" ~ "Quercus faginea/pyrenaica type",
                       variable == "Quercus pyrenaica" ~ "Quercus faginea/pyrenaica type",
                       variable == "Quercus ilex" ~ "Quercus ilex/coccifera type",
                       variable == "Quercus rotundifolia" ~ "Quercus ilex/rotundifolia/coccifera type",
                       variable == "Quercus coccifera" ~ "Quercus ilex/rotundifolia/coccifera type")) %>% 
    group_by(pixelId, variable, age) %>% 
    summarize(value = max(value)) %>% 
    ungroup()
  
  wide_neotoma <- neotoma_data %>% 
    pivot_wider(names_from = variablename,
                values_from = value,
                values_fill = 0)
  
  
  wide_pred <- pred_points %>% 
    pivot_wider(names_from = variable,
                values_from = value, 
                values_fill = 0)
  
  wide_neotoma <- wide_neotoma %>% 
    mutate(place_time = paste0(pixelId, "_", timebins)) %>% 
    mutate(period_block = cut(timebins,
                              breaks = c(-Inf, seq(1000, 21000, 1000), Inf),
                              labels = c(1:22),
                              right = TRUE))
  
  wide_pred <- wide_pred %>% 
    mutate(place_time = paste0(pixelId, "_", age)) %>% 
    mutate(period_block = cut(age,
                              breaks = c(-Inf, seq(1000, 21000, 1000), Inf),
                              labels = c(1:22),
                              right = TRUE))
  
  eval_results <- data.frame()
  
  periods <- wide_neotoma %>% pull(period_block) %>% unique() %>% as.numeric() %>% sort()
  
  for(j in 1:length(pollen_types)){
    
    for(i in 1:length(periods)){
      
      temp_neotoma <- wide_neotoma %>%
        select(place_time, period_block, pollen_types[j]) %>% 
        filter(period_block == periods[i])
      colnames(temp_neotoma) <- c("place_time", "period_block", "value")
      
      temp_pred <- wide_pred %>% 
        select(place_time, period_block, pollen_types[j]) %>% 
        filter(period_block == periods[i])
      colnames(temp_pred) <- c("place_time", "period_block", "value")
      
      
      p <- temp_neotoma %>% dplyr::filter(value == 1) %>% pull(place_time)
      a <- temp_neotoma %>% dplyr::filter(value == 0) %>% pull(place_time)

      index <- nrow(eval_results)
      
      if(length(p) < 1){
          print("no presences or low records")
        
        eval_results[index + 1, ] <- NA
        eval_results[index + 1, "pollen_types"] <- pollen_types[j]
        eval_results[index + 1, "period_block"] <- periods[i]
        
      }else{
        
        eval <- dismo::evaluate(temp_pred %>% dplyr::filter(place_time %in% p) %>% pull(value),
                                temp_pred %>% dplyr::filter(place_time %in% a) %>% pull(value))
        
        
        t <- dismo::threshold(eval)
        
        eval_results[index + 1, "AUC"] <- eval@auc
        eval_results[index + 1, "KAPPA"] <- eval@kappa[first(which(eval@t >= t$spec_sens))]
        eval_results[index + 1, "TSS"] <- eval@TPR[first(which(eval@t >= t$spec_sens))] + eval@TNR[first(which(eval@t >= t$spec_sens))] - 1
        eval_results[index + 1, "pollen_types"] <- pollen_types[j]
        eval_results[index + 1, "period_block"] <- periods[i]
      }
      
    }}
  
  eval_results
}


past.projection <- function(climate, soil, species, model, CM, period, i) {
  
  climate <- do.call(Map, c(c, climate))
  bio <- climate[[1]][[period]] %>% terra::unwrap() %>% subset("bio14")
  env <- climate[[2]][[period]] %>% terra::unwrap() %>% subset(c("aridityIndexThornthwaite", "continentality", "PETseasonality", "thermicityIndex"))
  names(env) <- c("aridityIndex", "continentality", "petSeasonality", "thermicityIndex")
  
  climate <- c(bio, env)
  soil <- soil %>% terra::unwrap()
  
  env_data <- c(climate, soil)
  
  if(CM == TRUE){
    XData.grid <- env_data %>% terra::values(dataframe = TRUE)
    xy.grid <- env_data %>% xyFromCell(cell = 1:ncell(.)) %>% as.data.frame()
    
    complete.index <- which(complete.cases(XData.grid))
    Gradient <- prepareGradient(model, XDataNew = XData.grid[complete.index,], sDataNew = list(units=xy.grid[complete.index,]))
    
    
    predY <- predict(model, Gradient = Gradient, expected = TRUE, predictEtaMean = TRUE, nParallel=2)
    
    EpredY <- matrix(NA, nrow(XData.grid), length(species))
    
    library(abind)
    EpredY[complete.index,] <- apply(abind(predY,along=3),c(1,2),mean)
    
    mapData <- data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY)
    colnames(mapData) <- c("x", "y", species)
    library(reshape2)
    mapData <- melt(mapData, id.vars = c("x", "y"))
    mapData$i <- i
    mapData$time <- period
  }else{
    
    pr <- sdm_predict(
      models = model,
      pred = env_data,
      thr = "max_sens_spec",
      con_thr = TRUE,
      predict_area = NULL
    )
    
    unconstrained <- pr$meanw[[1]]
    names(unconstrained) <- "value"
    
    mapData <- unconstrained %>% as.data.frame(xy = TRUE)
    mapData$i <- i
    mapData$variable <- species
    mapData$time <- period
    
  }
  
  mapData
  
}


weight.rast <- function(r, th) {
  r <- r %>% unwrap()
  
  recl_r <- ifel(r >= th, 1, NA)
  
  recl_r <- ifel(is.na(recl_r), 0, 1)
  
  
  mean_neighbors <- function(cell_values) {
    
    mean(cell_values[-221], na.rm = TRUE)  
  }
  
  
  focal_raster <- focal(recl_r, w = matrix(1, nrow = 21, ncol = 21), fun = mean_neighbors)
  
  weight_r <- 1 - focal_raster
  
  weight_r <- ifel(recl_r == 0, NA, weight_r)
  names(weight_r) <- "weight"
  names(r) <- "complection"
  
  output <- c(r, weight_r) %>% as.data.frame(cells = TRUE)
  colnames(output) <- c("Cell_id", "complection", "weight")
  output
}

slow.calibrate.hmsc <- function(data, thin, samples, transient, nChains, sp_names, var_names, coord_names, formula, nParallel){
  
  set.model <- function(data, sp_names, var_names, coord_names, formula) {
    Y <- data[, sp_names] 
    XData <- data[, var_names]
    xy <- data[, coord_names]
    Y <- as.matrix(Y)
    rownames(Y) <- rownames(XData)
    XFormula <- formula
    studyDesign <- data.frame(units = rownames(Y))
    studyDesign$units <- factor(studyDesign$units, levels = rownames(Y))
    rl <- HmscRandomLevel(sData = xy, sMethod = 'NNGP', nNeighbours = 10)
    rl = setPriors(rl, nfMin = 1, nfMax = 3)
    Hmsc(Y = Y,
         XData = XData,
         XFormula = XFormula,
         studyDesign = studyDesign,
         ranLevels = list("units" = rl),
         distr = "probit")
  }
  
  model.pa <- set.model(data, sp_names, var_names, coord_names, formula)
  
  # fitting the model
  model <- sampleMcmc(model.pa, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = nParallel)
  
  model
  
}

single.hmsc.calibration <- function(sp_names, data, thin, samples, transient, nchains, var_names, coord_names, formula, nParallel) {
  model_list <- list()
  
  for(j in 1:length(sp_names)){
    
    model <- slow.calibrate.hmsc(data, thin, samples, transient,
                            nchains, sp_names[j], var_names,
                            coord_names, formula, nParallel)
    
    model_list[[j]] <- model
    
  }
  
  model_list
}

single.past.projection <- function(sp, env, soil, model, CM, period, i) {
  projection_list <- list()
  
  for(z in 1:length(sp)){
    
    projection <- past.projection(env, soil, sp[z], model[[z]], CM = CM, period = period, i = i)  
    projection_list[[z]] <- projection
  }
  
  proj <- do.call(bind_rows, projection_list) 
  
  proj
}

eval.pred <- function(pred, obse, sp_names, clim_data, i, dataset) {
  
  pred <- pred  %>% 
    group_by(x, y, variable, time) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>%
    st_as_sf(coords = c("x", "y"))
  
  clim_data <- clim_data %>% terra::unwrap()
  
  template <- clim_data %>% subset(1)
  
  pred <- pred %>% mutate(pixelId = cellFromXY(template, pred %>% st_coordinates())) 
  
  coords <- template %>% xyFromCell(cell = pred$pixelId) %>% as.data.frame()
  
  pred <- bind_cols(pred, coords) %>% 
    st_drop_geometry()
  
  points <- pred %>% 
    dplyr::select(pixelId) |>  
    distinct(pixelId, .keep_all = TRUE) %>% 
    pull(pixelId)
  
  
  obse <- obse %>% 
    select(decimalLongitude, decimalLatitude, species) %>% 
    mutate(value = 1) %>% 
    pivot_wider(names_from = species,
                values_from = value, 
                values_fill = 0) %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"))
  
  obse <- obse %>% mutate(pixelId = cellFromXY(template, obse %>% st_coordinates())) 
  
  coords <- template %>% xyFromCell(cell = obse$pixelId) %>% as.data.frame()
  
  obse <- bind_cols(obse, coords) %>% 
    st_drop_geometry() %>% 
    filter(pixelId %in% points) %>% 
    group_by(pixelId, x, y) %>% 
    summarise_all("sum")
  
  obse[,sp_names] <- sapply(obse[,sp_names],function(x) ifelse(x >= 1, 1, 0))
  
  pred <- pred %>% filter(pixelId %in% obse$pixelId) %>% select(-pixelId)
  
  obse <- obse %>% ungroup() %>%  select(-pixelId)
  
  pred <- pred %>% 
    pivot_wider(names_from = variable,
                values_from = value)
  
  data <- obse %>% left_join(pred, by = c("x", "y"))
  
  thres <- data.frame(kappa = NA,
                      spec_sens = NA,
                      no_omission = NA,
                      prevalence = NA,
                      equal_sens_spec = NA,
                      sensitivity = NA,
                      AUC = NA,
                      KAPPA = NA,
                      TSS = NA)
  
  for(j in 1:length(sp_names)){
    print(j)
    
    eval <- dismo::evaluate(data %>% filter(get(paste0(sp_names[j], ".x"), pos = data) == 1) %>% pull(paste0(sp_names[j], ".y")),
                            data %>% filter(get(paste0(sp_names[j], ".x"), pos = data) == 0) %>% pull(paste0(sp_names[j], ".y")))
    thres[j,] <- dismo::threshold(eval)
    thres[j, "AUC"] <- eval@auc
    thres[j, "KAPPA"] <- eval@kappa[first(which(eval@t >= thres$spec_sens))]
    thres[j, "TSS"] <- eval@TPR[first(which(eval@t >= thres$spec_sens))] + eval@TNR[first(which(eval@t >= thres$spec_sens))] - 1
  }
  
  thres$sp <- sp_names
  
  thres
  
}

th.pred <- function(th, pred, sp_names) {
  
  th <- th %>%
    select(spec_sens, sp) 
  
 
  temp_list <- list()
  
  for(i in 1:length(sp_names)){
    
    th_value <- th %>%
      filter(sp == sp_names[i]) %>% 
      pull(spec_sens)
    
    temp_list[[i]] <- pred %>% 
      filter(variable == sp_names[i]) %>% 
      mutate(th_pred = as.numeric(value >= th_value))
  }
  
  names(temp_list) <- sp_names
  
  recla_pred <- do.call(bind_rows, temp_list)
  
  
  recla_pred
}

past.map<- function(clm_pred, sdm_pred, obse, sp_names, clim_data, name, th){
  
  clm_pred <- clm_pred %>% 
    group_by(x, y, variable, time) %>% 
    summarize(value = mean(value)) %>% 
    mutate(model = "CLM")
  
  sdm_pred <- sdm_pred %>% 
    group_by(x, y, variable, time) %>% 
    summarize(value = mean(value)) %>% 
    mutate(model = "SDM")
  
  recla_clm <- th.pred(th[[1]], clm_pred, sp_names) 
  recla_sdm <- th.pred(th[[2]], sdm_pred, sp_names)  
  
  recla_data <-  bind_rows(recla_clm, recla_sdm)
  
  clm <- recla_data %>% ungroup() %>% filter(model == "CLM") %>% mutate(th_pred = th_pred*10)
  recla_data <- recla_data %>% ungroup() %>%  filter(model == "SDM") %>% mutate(th_pred = th_pred + clm$th_pred)
  
   
  recla_data %>% ggplot() +
    geom_tile(aes(x = x, y = y, fill = as.factor(th_pred), col = as.factor(th_pred))) +
    scale_fill_manual(values = c("lightgrey", "yellow", "red", "orange" ), na.value="royalblue", 
                      labels = c("None", "SDM", "CLM", "Both")) +
    scale_colour_manual(values = c("lightgrey", "yellow", "red", "orange" ), na.value="royalblue") +
    guides(colour = "none") +
    facet_wrap(.~variable) +
    coord_equal() +
    theme_bw() +
    labs(fill = NULL)
  
  
  ggsave(paste0("./figures/", name, ".jpg"), height = 6, width = 12)
  
  recla_data
}

single.timeseries.prediction <- function(comm, vars, coords, species, model, i, dataset) {
  
  clim_list <- comm[[1]]
  soil_list <- comm[[2]]
  coords_list <- comm[[3]]
  
  data <- list()
  
  for(f in 1:length(clim_list)){
    data[[f]] <- bind_cols(clim_list[[f]], soil_list[[f]], coords_list[[f]])
  }
  
  pred_list <- list()
  
  for(z in 1:length(species)){
    sp <- species[[z]]
    mod <- model[[z]]
    
    prediction <- data |> lapply(FUN = function(x, y, z, a, b, c, w){
      coords.prediction(x, vars, coords, sp, mod, i, dataset)
    })
    
    prediction <- bind_rows(prediction, .id = "time")
    
    pred_list[[z]] <- prediction
    
  }
  
  output <- do.call(bind_rows, pred_list)
  output
  
}



##########


median2 = function(x) {
  return(median(x, na.rm = TRUE))
}

mean2 = function(x) {
  return(mean(x, na.rm = TRUE))
}

#Required for sensitivity analysis in performance metrics calculation

mean_prediction <- function(predY, hM) {
  
  mPredY = matrix(NA, nrow = nrow(predY), ncol = ncol(predY))
  sel = hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    mPredY[, sel] = as.matrix(apply(abind(predY[, sel, , 
                                                drop = FALSE], along = 3), c(1, 2), median2))
  }
  sel = !hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    mPredY[, sel] = as.matrix(apply(abind(predY[, sel, , 
                                                drop = FALSE], along = 3), c(1, 2), mean2))
  }
  colnames(mPredY) <- colnames(predY)
  rownames(mPredY) <- rownames(predY)
  mPredY
}



#modification of Hmsc function, making it more flexible for validation with external data

custom_evaluateModelFit <- function (obsY, predY, hM) {
  computeRMSE = function(obs, pred) {
    ns = dim(obs)[2]
    RMSE = rep(NA, ns)
    for (i in 1:ns) {
      RMSE[i] = sqrt(mean((obs[, i] - pred[, i])^2, na.rm = TRUE))
    }
    return(RMSE)
  }
  computeR2 = function(obs, pred, method = "pearson") {
    ns = dim(obs)[2]
    R2 = rep(NA, ns)
    for (i in 1:ns) {
      co = cor(obs[, i], pred[, i], method = method, use = "pairwise")
      R2[i] = sign(co) * co^2
    }
    return(R2)
  }
  computeAUC = function(obs, pred) {
    ns = dim(obs)[2]
    AUC = rep(NA, ns)
    obsY <- ifelse(obs > 0, 1, 0)
    for (i in 1:ns) {
      sel = !is.na(obs[, i])
      if (length(unique(obs[sel, i])) == 2) 
        AUC[i] = pROC::auc(obs[sel, i], pred[sel, i], levels = c(0, 
                                                                 1), direction = "<")
    }
    return(AUC)
  }
  computeTjurR2 = function(obs, pred) {
    ns = dim(obs)[2]
    R2 = rep(NA, ns)
    for (i in 1:ns) {
      R2[i] = mean(pred[which(obs[, i] == 1), i]) - mean(pred[which(obs[, 
                                                                        i] == 0), i])
    }
    return(R2)
  }
  
  mPredY = matrix(NA, nrow = nrow(predY), ncol = ncol(predY))
  sel = hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    mPredY[, sel] = as.matrix(apply(abind(predY[, sel, , 
                                                drop = FALSE], along = 3), c(1, 2), median2))
  }
  sel = !hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    mPredY[, sel] = as.matrix(apply(abind(predY[, sel, , 
                                                drop = FALSE], along = 3), c(1, 2), mean2))
  }
  RMSE = computeRMSE(obsY, mPredY)
  R2 = NULL
  AUC = NULL
  TjurR2 = NULL
  SR2 = NULL
  O.AUC = NULL
  O.TjurR2 = NULL
  O.RMSE = NULL
  C.SR2 = NULL
  C.RMSE = NULL
  sel = hM$distr[, 1] == 1
  if (sum(sel) > 0) {
    R2 = rep(NA, hM$ns)
    R2[sel] = computeR2(obsY[, sel, drop = FALSE], mPredY[, 
                                                          sel, drop = FALSE])
  }
  sel = hM$distr[, 1] == 2
  if (sum(sel) > 0) {
    AUC = rep(NA, hM$ns)
    TjurR2 = rep(NA, hM$ns)
    AUC[sel] = computeAUC(obsY[, sel, drop = FALSE], mPredY[, 
                                                            sel, drop = FALSE])
    TjurR2[sel] = computeTjurR2(obsY[, sel, drop = FALSE], 
                                mPredY[, sel, drop = FALSE])
  }
  sel = hM$distr[, 1] == 3
  if (sum(sel) > 0) {
    SR2 = rep(NA, hM$ns)
    O.AUC = rep(NA, hM$ns)
    O.TjurR2 = rep(NA, hM$ns)
    O.RMSE = rep(NA, hM$ns)
    C.SR2 = rep(NA, hM$ns)
    C.RMSE = rep(NA, hM$ns)
    SR2[sel] = computeR2(obsY[, sel, drop = FALSE], mPredY[, 
                                                           sel, drop = FALSE], method = "spearman")
    predO = 1 * (predY[, sel, , drop = FALSE] > 0)
    mPredO = as.matrix(apply(abind(predO, along = 3), c(1, 
                                                        2), mean2))
    O.AUC[sel] = computeAUC(1 * (obsY[, sel, drop = FALSE] > 
                                   0), mPredO)
    O.TjurR2[sel] = computeTjurR2(1 * (obsY[, sel, drop = FALSE] > 
                                         0), mPredO)
    O.RMSE[sel] = computeRMSE(1 * (obsY[, sel, drop = FALSE] > 
                                     0), mPredO)
    mPredCY = mPredY[, sel, drop = FALSE]/mPredO
    CY = obsY[, sel, drop = FALSE]
    CY[CY == 0] = NA
    C.SR2[sel] = computeR2(CY, mPredCY, method = "spearman")
    C.RMSE[sel] = computeRMSE(CY, mPredCY)
  }
  MF = list(RMSE = RMSE)
  if (!is.null(R2)) {
    MF$R2 = R2
  }
  if (!is.null(AUC)) {
    MF$AUC = AUC
  }
  if (!is.null(TjurR2)) {
    MF$TjurR2 = TjurR2
  }
  if (!is.null(SR2)) {
    MF$SR2 = SR2
  }
  if (!is.null(O.AUC)) {
    MF$O.AUC = O.AUC
  }
  if (!is.null(O.TjurR2)) {
    MF$O.TjurR2 = O.TjurR2
  }
  if (!is.null(O.RMSE)) {
    MF$O.RMSE = O.RMSE
  }
  if (!is.null(C.SR2)) {
    MF$C.SR2 = C.SR2
  }
  if (!is.null(C.RMSE)) {
    MF$C.RMSE = C.RMSE
  }
  return(MF)
}


#modification of Hmsc function, making it more flexible for validation with external data

tuning_computePredictedValues <- function(newData, oldData, modelo, Yc=NULL, vars) {
  require(dplyr)
  new_index <- !rownames(newData) %in% rownames(oldData)
  old_index <- rownames(newData) %in% rownames(oldData)
  
  old_sites <- rownames(newData)[old_index]
  rL1 = modelo$ranLevels[[1]]
  
  if(sum(old_index) > 0) {
    StudyDesignOld <- filter(modelo$studyDesign, units %in% old_sites)
    xyOld = filter(rL1$s, rownames(rL1$s) %in% old_sites)
  }
  
  if(sum(new_index) > 0) {
    com <- newData[new_index,]
    grad <- prepareGradient(modelo, XDataNew = com[, vars], sDataNew = list(units = com[, c("x", "y")]))
    xyNew = com[,1:2] 
    rownames(xyNew) = grad$StudyDesignNew[,1]
    colnames(xyNew) = colnames(xyOld)
  } 
  
  if(sum(new_index) > 0 & sum(old_index) > 0){
    StudyDesignAll <- rbind(StudyDesignOld, grad$studyDesignNew)
    xy <- rbind(xyOld, xyNew)
  }
  if(sum(new_index) > 0 & sum(old_index) == 0){
    StudyDesignAll <- grad$studyDesignNew
    xy <- xyNew
  }
  if(sum(new_index) == 0 & sum(old_index) > 0){
    StudyDesignAll <- StudyDesignOld
    xy <- xyOld
  }
  
  rL1$pi = StudyDesignAll[,1]
  rL1$s = xy
  rL1$N = nrow(rL1$s)
  
  if(is.null(Yc)){
    predYR1 = predict(modelo, studyDesign = StudyDesignAll, XData = newData[, vars], ranLevels = list("units"=rL1),
                      expected = TRUE, predictEtaMean = TRUE)
  }else{ 
    Yc[,1] <- NA
    predYR1 = predict.Hmsc(modelo, studyDesign = StudyDesignAll, XData = newData[, vars], ranLevels = list("units"=rL1),
                           expected = TRUE, predictEtaMean = TRUE, Yc=as.matrix(Yc))
    
  }
  predArray = abind(predYR1, along = 3)
  predArray
}


#sensitivity specificty analysis for testing model performance

sensibilidad_especifidad <- function(sp_names, pred, obs, tr) {
  
  pred <- as.data.frame(pred)
  pred$joinID <- rownames(pred)
  obs$joinID <- rownames(obs)
  
  original_cut<- semi_join(obs,pred, by= c("joinID"="joinID"))
  
  presencias_originalsp <- original_cut[which(original_cut[,sp_names]==1),]
  ausencias_originalsp <- original_cut[which(original_cut[,sp_names]==0),]
  
  presenciaspredsp<- pred
  presenciaspredsp<- presenciaspredsp[presencias_originalsp$joinID, ]
  presenciaspredsp<- presenciaspredsp[,sp_names] 
  
  ausenciaspredsp <- pred
  ausenciaspredsp <- ausenciaspredsp[ausencias_originalsp$joinID, ]
  ausenciaspredsp <- ausenciaspredsp[,sp_names] 
  
  e <- dismo::evaluate(p=presenciaspredsp, a=ausenciaspredsp, tr=tr)
  t <- threshold(e, "spec_sens")
  e <- dismo::evaluate(p=presenciaspredsp, a=ausenciaspredsp, tr=t)
  return(e)
}

extraccion <- function(e){
  e <- data.frame(e@TPR, e@TNR, e@kappa, e@t)
  
  colnames(e) <- c("TPR", "TNR", "kappa", "t")
  return(e)
}


tuning.metrics<- function(com_testing, com_training, model, sp_names, Yc=NULL, i, vars) {
  
  predY<- tuning_computePredictedValues(com_testing, com_training, model, Yc=Yc, vars)
  mpredY<- mean_prediction(predY, model)
  colnames(mpredY) <- sp_names
  metricas<- data.frame(custom_evaluateModelFit(obsY = data.frame(com_testing[, sp_names]), predY = predY, hM = model), species = sp_names)
  
  sensitivity <- lapply(sp_names, FUN= sensibilidad_especifidad, mpredY, com_testing)
  sensitivity <- do.call("rbind", lapply(sensitivity, FUN= extraccion))
  
  metricas <- cbind(metricas, sensitivity)
  metricas[,"TSS"]<- metricas$TPR + metricas$TNR - 1
  metricas[,"i"] <- i
  metricas <- reshape2::melt(metricas, measure.vars = c("RMSE", "AUC", "TjurR2", "TPR", "TNR", "kappa", "TSS"))
  metricas
}


#Gives performance metrics (threshold dependent and independent) for predictive power

tuning.predictive <- function(com_testing, model, sp_names, i) {
  predY = computePredictedValues(hM= model, partition= com_testing$partition , nParallel = 5)
  mpredY<- data.frame(mean_prediction(predY, model))
  colnames(mpredY)<- sp_names
  rownames(mpredY)<- rownames(com_testing)
  
  PP <- data.frame(evaluateModelFit(hM=model, predY=predY), species=sp_names)
  
  sensitivity <- lapply(sp_names, FUN= sensibilidad_especifidad, mpredY, com_testing)
  sensitivity <- do.call("rbind", lapply(sensitivity, FUN= extraccion))
  
  PP <- cbind(PP, sensitivity)
  PP[,"TSS"]<- PP$TPR + PP$TNR - 1
  PP[,"i"] <- i
  PP <- reshape2::melt(PP, measure.vars = c("RMSE", "AUC", "TjurR2", "TPR", "TNR", "kappa", "TSS"))
  PP
}

single.tuning.metrics <- function(sp_names, single_model, input, i_target, var_names) {
  exp_list <- list()
  
  for(z in 1:length(sp_names)){
    
    
    sp <- sp_names[[z]]
    mod <- single_model[[z]]
    
    exp <- tuning.metrics(input, input, mod, sp, Yc=NULL, i_target, var_names)
    
    exp_list[[z]] <- exp
    
  }
  
  output <- do.call(bind_rows, exp_list)
  output
}


single.tuning.predictive <- function(sp_names, single_model, input, i) {
  pred_list <- list()
  
  for(z in 1:length(sp_names)){
    
    
    sp <- sp_names[[z]]
    mod <- single_model[[z]]
    
    pred <- tuning.predictive(input, mod, sp, i) 
    
    
    pred_list[[z]] <- pred
  }
  
  output <- do.call(bind_rows, pred_list)
  output
}

single.spatial.prediction <- function(climate, soil, species, model, i, dataset) {
  
  prediction_list <- list()
  
  for(z in 1:length(species)){
    
    prediction <- spatial.prediction(climate, soil, species[z], model[[z]], i, dataset)  
    prediction_list[[z]] <- prediction
  }
  
  pred <- do.call(bind_rows, prediction_list) 
  
  pred
}

eval.pred.pres <- function(pred, obse, sp_names, clim_data, i, dataset) {
  
  pred <- pred  %>% 
    st_as_sf(coords = c("x", "y"))
  
  clim_data <- clim_data %>% terra::unwrap()
  
  template <- clim_data %>% subset(1)
  
  pred <- pred %>% mutate(pixelId = cellFromXY(template, pred %>% st_coordinates())) 
  
  coords <- template %>% xyFromCell(cell = pred$pixelId) %>% as.data.frame()
  
  pred <- bind_cols(pred, coords) %>% 
    st_drop_geometry()
  
  points <- pred %>% 
    dplyr::select(pixelId) |>  
    distinct(pixelId, .keep_all = TRUE) %>% 
    pull(pixelId)
  
  
  obse <- obse %>% 
    select(decimalLongitude, decimalLatitude, species) %>% 
    mutate(value = 1) %>% 
    pivot_wider(names_from = species,
                values_from = value, 
                values_fill = 0) %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"))
  
  obse <- obse %>% mutate(pixelId = cellFromXY(template, obse %>% st_coordinates())) 
  
  coords <- template %>% xyFromCell(cell = obse$pixelId) %>% as.data.frame()
  
  obse <- bind_cols(obse, coords) %>% 
    st_drop_geometry() %>% 
    filter(pixelId %in% points) %>% 
    group_by(pixelId, x, y) %>% 
    summarise_all("sum")
  
  obse[,sp_names] <- sapply(obse[,sp_names],function(x) ifelse(x >= 1, 1, 0))
  
  pred <- pred %>% filter(pixelId %in% obse$pixelId) %>% select(-pixelId)
  
  obse <- obse %>% ungroup() %>%  select(-pixelId)
  
  pred <- pred %>% 
    pivot_wider(names_from = variable,
                values_from = value)
  
  data <- obse %>% left_join(pred, by = c("x", "y"))
  
  thres <- data.frame(kappa = NA,
                      spec_sens = NA,
                      no_omission = NA,
                      prevalence = NA,
                      equal_sens_spec = NA,
                      sensitivity = NA,
                      AUC = NA,
                      KAPPA = NA,
                      TSS = NA)
  
  for(j in 1:length(sp_names)){
    print(j)
    
    eval <- dismo::evaluate(data %>% filter(get(paste0(sp_names[j], ".x"), pos = data) == 1) %>% pull(paste0(sp_names[j], ".y")),
                            data %>% filter(get(paste0(sp_names[j], ".x"), pos = data) == 0) %>% pull(paste0(sp_names[j], ".y")))
    thres[j,] <- dismo::threshold(eval)
    thres[j, "AUC"] <- eval@auc
    thres[j, "KAPPA"] <- eval@kappa[first(which(eval@t >= thres$spec_sens))]
    thres[j, "TSS"] <- eval@TPR[first(which(eval@t >= thres$spec_sens))] + eval@TNR[first(which(eval@t >= thres$spec_sens))] - 1
  }
  
  thres$sp <- sp_names
  
  thres
  
}


get.threshold <- function(clm_pred, sdm_pred, obse, sp_names, clim_data) {
  clm_pred <- clm_pred %>% 
    group_by(x, y, variable) %>% 
    summarize(value = mean(value)) %>% 
    mutate(model = "CLM")
  
  sdm_pred <- sdm_pred %>% 
    group_by(x, y, variable) %>% 
    summarize(value = mean(value)) %>% 
    mutate(model = "SDM")
  
  
  clm_th <- eval.pred.pres(clm_pred, obse, sp_names, clim_data) 
  sdm_th <- eval.pred.pres(sdm_pred, obse, sp_names, clim_data) 
  
  
  th <- list()
  
  th[[1]] <- clm_th
  th[[2]] <- sdm_th
  
  th
}


present.map<- function(clm_pred, sdm_pred, obse, sp_names, clim_data, name, th){
  
  clm_pred <- clm_pred %>% 
    group_by(x, y, variable) %>% 
    summarize(value = mean(value)) %>% 
    mutate(model = "CLM")
  
  sdm_pred <- sdm_pred %>% 
    group_by(x, y, variable) %>% 
    summarize(value = mean(value)) %>% 
    mutate(model = "SDM")
  
  recla_clm <- th.pred(th[[1]], clm_pred, sp_names) 
  recla_sdm <- th.pred(th[[2]], sdm_pred, sp_names)  
  
  recla_data <-  bind_rows(recla_clm, recla_sdm)  
  
  clm <- recla_data %>% ungroup() %>% filter(model == "CLM") %>% mutate(th_pred = th_pred*10)
  recla_data <- recla_data %>% ungroup() %>%  filter(model == "SDM") %>% mutate(th_pred = th_pred + clm$th_pred)
  
  
  recla_data %>% ggplot() +
    geom_tile(aes(x = x, y = y, fill = as.factor(th_pred), col = as.factor(th_pred))) +
    scale_fill_manual(values = c("lightgrey", "yellow", "red", "orange" ), na.value="royalblue", 
                      labels = c("None", "SDM", "CLM", "Both")) +
    scale_colour_manual(values = c("lightgrey", "yellow", "red", "orange" ), na.value="royalblue") +
    guides(colour = "none") +
    facet_wrap(.~variable) +
    coord_equal() +
    theme_bw() +
    labs(fill = NULL)
  
  ggsave(paste0("./figures/", name, ".jpg"), height = 6, width = 12)
  
  recla_data
}


refugia.map <- function(present, MH, LGM, name){
  recla <- function(x, model = "SDM") {
    
    if(model == "CLM"){
      x <- x %>% 
        mutate(th_pred = case_when(
          th_pred == 10 ~ 1,
          th_pred == 1 ~ 0, 
          th_pred == 11 ~ 1, 
          .default = th_pred
        ))
    }else{
      x <- x %>% 
        mutate(th_pred = case_when(
          th_pred == 10 ~ 0,
          th_pred == 1 ~ 1, 
          th_pred == 11 ~ 1, 
          .default = th_pred
        ))
    }
    x  
    
  }
  
  
  CLM_present <- present %>% recla("CLM")
  CLM_MH <- MH %>% recla("CLM")
  CLM_LGM <- LGM %>% recla("CLM")
  
  CLM <- CLM_present %>% 
    mutate(th_pred = th_pred + CLM_MH$th_pred + CLM_LGM$th_pred) %>% 
    mutate(th_pred = case_when(
      th_pred == 1 ~ 0, 
      th_pred == 2 ~ 0, 
      th_pred == 3 ~ 10, 
      .default = th_pred
    ))
  
  SDM_present <- present %>% recla()
  SDM_MH <- MH %>% recla()
  SDM_LGM <- LGM %>% recla()
  
  SDM <- SDM_present %>% 
    mutate(th_pred = th_pred + SDM_MH$th_pred + SDM_LGM$th_pred) %>% 
    mutate(th_pred = case_when(
      th_pred == 1 ~ 0, 
      th_pred == 2 ~ 0, 
      th_pred == 3 ~ 1, 
      .default = th_pred
    ))
  
  data <- SDM %>% 
    mutate(th_pred = th_pred + CLM$th_pred)
  
  data %>% ggplot() +
    geom_tile(aes(x = x, y = y, fill = as.factor(th_pred), col = as.factor(th_pred))) +
    scale_fill_manual(values = c("lightgrey", "yellow", "red", "orange" ), na.value="royalblue", 
                      labels = c("None", "SDM", "CLM", "Both")) +
    scale_colour_manual(values = c("lightgrey", "yellow", "red", "orange" ), na.value="royalblue") +
    guides(colour = "none") +
    facet_wrap(.~variable) +
    coord_equal() +
    theme_bw() +
    labs(fill = NULL)
  
  ggsave(paste0("./figures/", name, ".jpg"), height = 6, width = 12)
  
}

prepare.future.data <- function(start, end, rcp, gcm){
  
  
  stars_to_stack <- function(i, data){
    data %>% slice("time", i:(i+11)) %>% 
      rast() %>% stack()
  }
  
  
  tmax <- dsclimtools::read_dsclim("./dsclim", "tasmax", start, end, calendar_dates = TRUE,
                                   proxy = FALSE, rcp = rcp, gcm = gcm)
  
  
  tmax <-  lapply(seq(1, (end-start+1)*12, by = 12), 
                  FUN = stars_to_stack, 
                  tmax) %>% 
    lapply(FUN = function(x){
      crs(x) <- CRS('+init=EPSG:4326')
      x
    })
  
  
  tmin <- dsclimtools::read_dsclim("./dsclim", "tasmin", start, end, calendar_dates = TRUE,
                                   proxy = FALSE, rcp = rcp, gcm = gcm)
  
  
  tmin <-  lapply(seq(1, (end-start+1)*12, by = 12), 
                  FUN = stars_to_stack, 
                  tmin) %>% 
    lapply(FUN = function(x){
      crs(x) <- CRS('+init=EPSG:4326')
      x
    })
  
  
  t <- dsclimtools::read_dsclim("./dsclim", "tas", start, end, calendar_dates = TRUE,
                                proxy = FALSE, rcp = rcp, gcm = gcm)
  
  
  t <-  lapply(seq(1, (end-start+1)*12, by = 12), 
               FUN = stars_to_stack, 
               t) %>% 
    lapply(FUN = function(x){
      crs(x) <- CRS('+init=EPSG:4326')
      x
    })
  
  
  prec <- dsclimtools::read_dsclim("./dsclim", "pr", start, end, calendar_dates = TRUE,
                                   proxy = FALSE, rcp = rcp, gcm = gcm)
  
  
  prec <-  lapply(seq(1, (end-start+1)*12, by = 12), 
                  FUN = stars_to_stack, 
                  prec) %>% 
    lapply(FUN = function(x){
      crs(x) <- CRS('+init=EPSG:4326')
      x
    })
  
  
  bio_list<- list()
  
  for(j in 1:length(prec)){
    
    bio_list[[j]] <- dismo::biovars(prec[[j]], tmin[[j]], tmax[[j]]) |> rast()
  }
  etsolrad_seq <- start:end
  
  etsolrad <- lapply(etsolrad_seq,
                     FUN = function(x, y){envirem::ETsolradRasters(y, year = x)},
                     prec[[1]] |> subset(1) |> rast())
  
  etsolrad <- do.call(mean, etsolrad)
  
  clim_list <- list()
  
  for(c in 1:length(prec)){
    
    temp_prec <- prec[[c]] |> rast()
    terra::set.names(temp_prec, paste0("precip_", 1:12))
    
    temp_t <- t[[c]] |> rast()
    terra::set.names(temp_t, paste0("tmean_", 1:12))
    
    temp_tmin <- tmin[[c]] |> rast()
    terra::set.names(temp_tmin, paste0("tmin_", 1:12))
    
    temp_tmax <- tmax[[c]] |> rast()
    terra::set.names(temp_tmax, paste0("tmax_", 1:12))
    
    clim_list[[c]] <- c(temp_prec, temp_t, temp_tmin, temp_tmax)
  }
  
  envirem_data <- lapply(clim_list, FUN = function(x){envirem::generateEnvirem(x, etsolrad, var = "all")}) |> 
    lapply(FUN = function(x){terra::wrap(x)})
  
  bio_list <- bio_list |> lapply(FUN = function(x){terra::wrap(x)})
  
  
  variables <- list()
  variables[[1]] <- bio_list
  variables[[2]] <- envirem_data
  variables
}

fut.projection <- function(fut_data, soil_data, species, model, period, i, rcp, gcm) {
  
  bio <- fut_data[[1]]
  env <- fut_data[[2]]
  
  bio_subset <- lapply(bio[period], FUN = unwrap) %>% do.call(mean, .) %>% subset("bio14")
  env_subset <- lapply(env[period], FUN = unwrap) %>% do.call(mean, .) %>% subset(c("aridityIndexThornthwaite", "continentality", "PETseasonality", "thermicityIndex"))
  names(env_subset) <- c("aridityIndex", "continentality", "petSeasonality", "thermicityIndex")
  
  
  fut_clim <-  c(bio_subset, env_subset)
  fut_clim %>% plot()
  
  soil_data <- soil_data %>% unwrap()
  
  env_data <- c(fut_clim, soil_data)
  
  
  XData.grid <- env_data %>% terra::values(dataframe = TRUE)
  xy.grid <- env_data %>% xyFromCell(cell = 1:ncell(.)) %>% as.data.frame()
  
  complete.index <- which(complete.cases(XData.grid))
  Gradient <- prepareGradient(model, XDataNew = XData.grid[complete.index,], sDataNew = list(units=xy.grid[complete.index,]))
  
  
  predY <- predict(model, Gradient = Gradient, expected = TRUE, predictEtaMean = TRUE, nParallel=2,
                   post = poolMcmcChains(model$postList)[seq(1, 200, by = 4)])
  
  EpredY <- matrix(NA, nrow(XData.grid), length(species))
  
  library(abind)
  EpredY[complete.index,] <- apply(abind(predY,along=3),c(1,2),mean)
  
  mapData <- data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY)
  colnames(mapData) <- c("x", "y", species)
  library(reshape2)
  mapData <- melt(mapData, id.vars = c("x", "y"))
  mapData$i <- i
  mapData$rcp <- rcp
  mapData$gcm <- gcm
  mapData$period <- paste0(first(period), "_", last(period))
  
  mapData
  
}

single.fut.projection <- function(fut_data, soil_data, species, model, period, i, rcp, gcm) {
  projection_list <- list()
  
  for(z in 1:length(species)){
    
    projection <- fut.projection(fut_data, soil_data, species[z], model[[z]], period, i, rcp, gcm)  
    
    projection_list[[z]] <- projection
  }
  
  proj <- do.call(bind_rows, projection_list) 
  
  proj
}

future.map <- function(fut_pred, name) {
  fut_pred <- fut_pred %>% 
    group_by(x, y, variable, rcp) %>% 
    summarize(value = mean(value))
  
  fut_pred %>% ggplot() +
    geom_tile(aes(x = x, y = y, col = value, fill = value)) +
    coord_equal() +
    facet_grid(variable~rcp)
  
  ggsave(paste0("./figures/", name, ".jpg"),  width = 6, height = 18)     
  
  fut_pred
}


eval.past <- function(CLM_eval, SDM_eval) {
  CLM <- CLM_eval %>%
    select(AUC, TSS, pollen_types, period_block) %>%
    pivot_longer(cols = c("AUC", "TSS"),
                 names_to = "metric",
                 values_to = "value"
    ) %>%
    mutate(model = "JSDM")
  
  SDM <- SDM_eval %>%
    select(AUC, TSS, pollen_types, period_block) %>%
    pivot_longer(cols = c("AUC", "TSS"),
                 names_to = "metric",
                 values_to = "value"
    ) %>%
    mutate(model = "SDM")
  
  data <- bind_rows(CLM, SDM)
  
  data <- data %>% 
    mutate(short_pollen_types = case_when(
      pollen_types == "Pinus halepensis/pinea type" ~ "P. halepens./pinea type",
      pollen_types == "Pinus nigra/sylvestris type" ~ "P. nigra/sylvestris. type",
      pollen_types == "Quercus faginea/pyrenaica type" ~ "Q. faginea./pyren. type",
      pollen_types == "Quercus ilex/rotundifolia/coccifera type" ~ "Q. ilex/rotun./coccif. type",
      .default = pollen_types
    ))
  
  
  library(ggpubr)
  
  
  data <- data %>% 
    mutate(metric = as.factor(metric))
  
  AUC1 <- levels(data$metric)[1]
  TSS1 <- levels(data$metric)[2]
  
  data <- data %>% 
    mutate(model = factor(model, levels = c("SDM", "JSDM")))
  
  data %>% ggplot() +
    geom_boxplot(aes(x = model, y = value, fill = model)) +
    facet_grid(metric~short_pollen_types) +
    facetted_pos_scales(
      y = list(
        metric == AUC1 ~ scale_y_continuous(limits = c(-0.2, 1)),
        metric == TSS1 ~ scale_y_continuous(limits = c(-0.2, 1))
      )
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 25),
          legend.title = element_blank(),
          legend.position = "bottom",
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          axis.text=element_text(size=15),
          legend.text = element_text(size = 25))  +
    stat_compare_means(aes(x = model, y = value, label = paste0("WSR. p = ", ..p.signif..)),
                       method = "wilcox.test", paired = TRUE, label.y.npc = 0, size = 6) +
    xlab("") +
    ylab("") 
  
  ggsave("./figures/past_eval_V4.jpg", height = 7, width = 16)
  
  data
}

present.eval <- function(CLM, SDM, name) {
  CLM <- CLM %>%
    filter(variable %in% c("AUC", "TSS")) %>%
    mutate(model = "JSDM")
  
  SDM <- SDM %>%
    filter(variable %in% c("AUC", "TSS")) %>%
    mutate(model = "SDM")
  
  data <- bind_rows(CLM, SDM)
  
  
  data %>% ggplot() +
    geom_boxplot(aes(x = model, y = value, fill = species)) +
    facet_grid(variable~species) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    stat_compare_means(aes(x = model, y = value), method = "t.test")
  
  ggsave(paste0("./figures/", name, ".jpg"), height = 7, width = 15)
  
  data
}

past.projection.loop <- function(past_env_data, soil, species, model, period, i) {
  
  
  proj_list <- list()
  
  for(p in 1:length(period)){
    climate <- do.call(Map, c(c, past_env_data))
    bio <- climate[[1]][[period[p]]] %>% terra::unwrap() %>% subset("bio14")
    env <- climate[[2]][[period[p]]] %>% terra::unwrap() %>% subset(c("aridityIndexThornthwaite", "continentality", "PETseasonality", "thermicityIndex"))
    names(env) <- c("aridityIndex", "continentality", "petSeasonality", "thermicityIndex")
    
    climate <- c(bio, env)
    soil <- soil %>% terra::unwrap()
    
    env_data <- c(climate, soil)
    
    XData.grid <- env_data %>% terra::values(dataframe = TRUE)
    xy.grid <- env_data %>% xyFromCell(cell = 1:ncell(.)) %>% as.data.frame()
    
    complete.index <- which(complete.cases(XData.grid))
    Gradient <- prepareGradient(model, XDataNew = XData.grid[complete.index,], sDataNew = list(units=xy.grid[complete.index,]))
    
    
    predY <- predict(model, Gradient = Gradient, expected = TRUE, predictEtaMean = TRUE, nParallel=2,
                     post = poolMcmcChains(model$postList)[seq(1, 200, by = 4)])
    
    
    EpredY <- matrix(NA, nrow(XData.grid), length(species))
    
    library(abind)
    EpredY[complete.index,] <- apply(abind(predY,along=3),c(1,2),mean)
    
    mapData <- data.frame(x=xy.grid[,1],y=xy.grid[,2], EpredY)
    colnames(mapData) <- c("x", "y", species)
    library(reshape2)
    mapData <- melt(mapData, id.vars = c("x", "y"))
    mapData$i <- i
    mapData$time <- period[p]
    
    proj_list[[p]] <- mapData
    
    rm(climate, bio, env, env_data, XData.grid, xy.grid, Gradient, predY, EpredY, mapData)
    gc()
    
  }
  
  proj_list
  
}

single.past.projection.loop <- function(env, soil, sp, model, period, i) {
  projection_list <- list()
  
  for(z in 1:length(sp)){
    
    projection <- past.projection.loop(env, soil, sp[z], model[[z]], period, i) 
    projection_list[[z]] <- projection
    
    gc()
  }
  
  proj <- do.call(bind_rows, projection_list) 
  
  proj
}

map.slope <- function(data, CM, th, myspecies, name) {
  data <- data %>% 
    group_by(x, y, variable, time) %>% 
    summarize(value = mean(value))
  
  
  slope <- function(x, t) {
    if (anyNA(x)){
      NA_real_
    } else {
      lm.fit(cbind(1, t), x)$coefficients[2]
    } 
  }
  
  
  result <- data %>%
    group_by(x, y, variable) %>%
    summarise(slope = slope(value, time), .groups = "drop")
  
  
  if(CM == TRUE){
    th <- th[[1]]
  }else{
    th <- th[[2]]
  }
  
  mask <- data
  
  mask <- th.pred(th, mask, myspecies) 
  
  mask <- mask %>% 
    group_by(x, y , variable) %>% 
    summarize(th_pred = max(th_pred))
  
  
  df_comb <- left_join(mask, result, by = c("x", "y", "variable"))
  
  df_comb <- df_comb %>%
    mutate(color_value = ifelse(th_pred == 1, slope, NA),
           is_gray = ifelse(th_pred == 0, TRUE, FALSE))
  
  ggplot(df_comb) +
    geom_tile(data = filter(df_comb, is_gray), aes(x = x, y = y), fill = "grey80") +
    geom_tile(data = filter(df_comb, !is_gray), aes(x = x, y = y, fill = color_value)) +
    facet_wrap(~ variable) +
    scale_fill_gradient2(
      name = "Value",
      low = "red",       # valores negativos
      mid = "lightyellow",    # cero
      high = "darkgreen",    # valores positivos
      midpoint = 0,
      na.value = "darkblue"  # NA de 'value'
    ) +
    coord_equal() +
    theme_minimal() 
  
  ggsave(paste0("./figures/", name, ".jpg"), width = 12, height = 6)
  
  df_comb
}

combine.eval <- function(exp_present, pred_present) {
  exp <- exp_present %>% 
    mutate(eval = "EXP") %>% 
    mutate(group = paste0(eval, "_", model))
  
  pred <- pred_present %>% 
    mutate(eval = "PRED") %>% 
    mutate(group = paste0(eval, "_", model))
  
  data <- bind_rows(exp, pred) %>% 
    mutate(short_species = case_when(
      species == "Abies alba" ~ "A. alba",
      species == "Abies pinsapo" ~ "A. pinsapo",
      species == "Cedrus atlantica" ~ "C. atlantica",
      species == "Pinus halepensis" ~ "P. halepens.",
      species == "Pinus nigra" ~ "P. nigra",
      species == "Pinus pinea" ~ "P. pinea",
      species == "Pinus sylvestris" ~ "P. sylvestris",
      species == "Quercus coccifera" ~ "Q. coccifera",
      species == "Quercus faginea" ~ "Q. faginea",
      species == "Quercus ilex" ~ "Q. ilex",
      species == "Quercus pyrenaica" ~ "Q. pyrenaica",
      species == "Quercus rotundifolia" ~ "Q. rotundif."
    ))
  
  
  exp <- data %>% 
    filter(eval == "EXP") %>% 
    droplevels() %>% 
    mutate(model = factor(model, levels = c("SDM", "JSDM")))
  
  pred <- data %>% 
    filter(eval == "PRED") %>% 
    droplevels() %>% 
    mutate(model = factor(model, levels = c("SDM", "JSDM")))
  
  AUC1 <- levels(exp$variable)[1]
  TSS1 <- levels(exp$variable)[2]
  AUC2 <- levels(pred$variable)[1]
  TSS2 <- levels(pred$variable)[2]
  
  
  plot1<- exp %>%
    ggplot(aes(x = model, y = value, fill = model)) +
    geom_boxplot() +
    facet_grid(rows = vars(variable), cols = vars(short_species), scales = "free_y") +
    facetted_pos_scales(
      y = list(
        variable == AUC1 ~ scale_y_continuous(limits = c(0.72, 1)),
        variable == TSS1 ~ scale_y_continuous(limits = c(0.35, 1))
      )
    ) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 17),
          legend.title = element_blank(),
          legend.position = "bottom",
          strip.text.x = element_text(size = 18),
          strip.text.y = element_text(size = 18),
          legend.text = element_text(size = 25)) +
    stat_compare_means(aes(x = model, y = value, label = paste0("WSR. p = ", ..p.signif..)), method = "wilcox.test", paired = TRUE,
                       label.y.npc = 0, size = 6) +
    xlab("") +
    ylab("")
  
  plot2<- pred %>%
    ggplot(aes(x = model, y = value, fill = model)) +
    geom_boxplot() +
    facet_grid(rows = vars(variable), cols = vars(short_species), scales = "free_y") +
    facetted_pos_scales(
      y = list(
        variable == AUC2 ~ scale_y_continuous(limits = c(0.72, 1)),
        variable == TSS2 ~ scale_y_continuous(limits = c(0.35, 1))
      )
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 25),
          axis.text.y = element_text(size = 17),
          legend.title = element_blank(),
          legend.position = "bottom",
          strip.text.x = element_text(size = 18),
          strip.text.y = element_text(size = 18),
          legend.text = element_text(size = 25)) +
    stat_compare_means(aes(x = model, y = value, label = paste0("WSR. p = ", ..p.signif..)), method = "wilcox.test", paired = TRUE,
                       label.y.npc = 0, size = 6) +
    xlab("") +
    ylab("")
  
  
  
  theme_set(theme_pubr())
  
  ggarrange(plot1, plot2, ncol = 1, nrow = 2, heights = c(6, 7.5),  align = "v",
            common.legend = TRUE, legend = "bottom", labels = "auto") %>%
    ggexport(filename = "./figures/combined_present_eval_v6.pdf", width = 20, height = 10)
  
  data
}

subset.future <- function(fut_map, name) {
  fut <- fut_map %>% 
    filter(variable %in% c("Cedrus atlantica", "Pinus nigra", 
                           "Quercus rotundifolia")) %>% 
    mutate(short_species = case_when(
      variable == "Cedrus atlantica" ~ "C. atlantica",
      variable == "Pinus nigra" ~ "P. nigra",
      variable == "Quercus rotundifolia" ~ "Q. rotundifolia"
    ))
  
  fut %>% ggplot() +
    geom_tile(aes(x = x, y = y, col = value, fill = value)) +
    scale_fill_viridis_c() +
    scale_colour_viridis_c() +
    theme_bw() +
    coord_equal() +
    facet_grid(short_species~rcp) +
    theme(legend.position = "bottom")
  
  ggsave(paste0("./figures/", name, ".jpg"), width = 12, height = 6)    
  
  fut
}


plot.arrange <- function(pollen_data, mask, all_timeseries_pred, complection_rast, indf) {
  neotoma.coords <- function(neotoma, mask, timeseries_pred) {
    pollen_types <- neotoma %>% colnames() %>% .[seq_along(.) %in% c(4:9)]
    
    neo_points <- neotoma  %>%  
      pivot_longer(cols = pollen_types,
                   names_to = "variablename",
                   values_to = "value") %>% 
      st_as_sf(coords = c("long", "lat"))
    
    neo_points <- neo_points %>% mutate(pixelId = cellFromXY(mask, neo_points %>% st_coordinates())) 
    
    neotoma_data <- neo_points %>% 
      st_drop_geometry() %>% 
      group_by(pixelId, variablename, timebins) %>% 
      summarize(value = max(value)) %>% 
      ungroup() 
    
    
    pred_points <- timeseries_pred   %>%  
      na.omit() %>%
      st_as_sf(coords = c("x", "y"))
    
    pred_points <- pred_points %>% mutate(pixelId = cellFromXY(mask, pred_points %>% st_coordinates())) %>% st_drop_geometry()
    
    pred_points <- pred_points %>% select(-dataset) %>% group_by(pixelId, variable, time) %>% 
      summarize(value = mean(value)) %>% ungroup() 
    
    pred_pixels <- pred_points %>% pull(pixelId) %>% unique()
    
    neotoma_pixels <- neotoma_data %>% pull(pixelId) %>% unique()
    
    pixel2keep <- intersect(pred_pixels, neotoma_pixels)
    
    neotoma_data <- neotoma_data %>% dplyr::filter(pixelId %in% pixel2keep) %>% as.data.frame()
    
    records_coords <- mask %>% xyFromCell(cell = neotoma_data$pixelId) %>% as.data.frame() %>% distinct()
    
    records_coords
  }
  
  
  neo_coords <- neotoma.coords(pollen_data, mask, all_timeseries_pred)
  
  
  nas <- mask %>% 
    as.data.frame(xy = TRUE, na.rm = FALSE)
  
  
  
  comp <- complection_rast %>%
    terra::unwrap() %>% 
    as.data.frame(xy = TRUE)
  
  colnames(comp) <- c("x", "y", "Completeness")
  
  
  
  #####below 10 records
  
  indf <- as.data.frame(indf)
  
  coords <- indf %>% 
    select(Longitude, Latitude) 
  
  indf$cust_cell_id <- cellFromXY(mask, coords)
  
  dat1 <- indf %>% 
    select(Scientific_name, Date_collected, cust_cell_id) %>% 
    distinct(Scientific_name, Date_collected, cust_cell_id)
  dat2 <- dat1 %>% group_by(cust_cell_id) %>% 
    summarize(cell_ct = n())
  dat3 <- dat2 %>% 
    filter(cell_ct < 10) 
  
  
  xy <- mask %>% 
    xyFromCell(cell = dat3$cust_cell_id) %>%
    as.data.frame() 
  
  
  pal <- c(
    "#f1b6da",  # rosa pastel
    "#ffffbf",  # amarillo crema (centro claro)
    "#e6f5d0",  # verde muy claro
    "#b8e186",  # verde medio
    "#4d9221"   # verde oscuro
  )
  
  comp_plot <- ggplot() +
    geom_tile(data = nas, aes(x = x, y = y, fill = bio1), show.legend = FALSE) +
    scale_fill_gradient(low = "lightgray", high = "lightgray", na.value = "darkblue") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = xy, aes(x = x, y = y), fill = pal[1]) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = comp, aes(x = x, y = y, fill = Completeness)) +
    scale_fill_gradientn(
      colours = pal[2:5],  
      name = "Comp."
    ) +
    theme_bw() +
    theme(
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.3)
    ) +
    labs(x = "Longitude", y = "Latitude") +
    coord_equal()
  
  
  records_plot <- ggplot() +
    geom_tile(data = nas, aes(x = x, y = y, fill = bio1), show.legend = FALSE) +
    scale_fill_gradient(low = "lightgray", high = "lightgray", na.value = "darkblue") +
    ggnewscale::new_scale_fill() +
    geom_point(data = neo_coords, aes(x = x, y = y, color = "Pollen records"), size = 2) +
    scale_color_manual(values = c("Pollen records" = "red")) +
    labs(color = "") +   
    theme_bw() +
    theme(legend.position = "none",
          legend.text = element_text(size = 18),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15)) +
    labs(x = "Longitude", y = "Latitude") +
    coord_equal() 
  
  
  sites_data <- pollen_data %>% 
    mutate(period_block = cut(timebins,
                              breaks = c(-Inf, seq(1000, 21000, 1000), Inf),
                              labels = c(1:22),
                              right = TRUE)) %>% 
    
    group_by(period_block, long, lat) %>% 
    summarize(Sites = 1)  %>% 
    group_by(period_block) %>% 
    summarize(Sites = sum(Sites)) %>% 
    ungroup() %>% 
    mutate(period_block = as.numeric(period_block))
  
  records_n_plot <- sites_data %>%
    ggplot() +
    geom_col(aes(x = period_block, y = Sites), fill = "darkblue", linewidth = 1.2) +
    xlab("Years (ka BP)") +
    ylab("Sites number") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15))
  
  
  bottom_row <- ggarrange(
    records_plot, records_n_plot,
    ncol = 2, labels = c("b", "c")
  )
  
  arr <- ggarrange(
    comp_plot,        
    bottom_row,       
    ncol = 1,
    heights = c(2, 1),   
    labels = c("a", "")) 
  
  ggsave(
    filename = "./figures/arrange_figure_V6.png",
    plot = arr,
    width = 8, height = 8, units = "in",
    dpi = 300,
    device = ragg::agg_png,
    bg = "white"
  )
}


plot.diff <- function(fut_subset, single_fut_subset, th, myspecies) {
  
  jsdm <- fut_subset %>% 
    mutate(jsdm = value) %>% 
    select(-value)
  
  sdm <- single_fut_subset %>% 
    mutate(sdm = value) %>% 
    select(-value)
  
  join <- full_join(jsdm, sdm) %>% 
    mutate("JSDM - SDM" = jsdm - sdm) %>% 
    mutate(rcp = case_when(
      rcp == "rcp2.6" ~ "RCP 2.6",
      rcp == "rcp4.5" ~ "RCP 4.5",
      rcp == "rcp6.0" ~ "RCP 6.0",
      rcp == "rcp8.5" ~ "RCP 8.5",
      .default = rcp
    ))
  
  
  recla_sdm <- th.pred(th[[2]], single_fut_subset, myspecies) %>% 
    mutate(rcp = case_when(
      rcp == "rcp2.6" ~ "RCP 2.6",
      rcp == "rcp4.5" ~ "RCP 4.5",
      rcp == "rcp6.0" ~ "RCP 6.0",
      rcp == "rcp8.5" ~ "RCP 8.5",
      .default = rcp
    ))
  
  
  
  ggplot() +
    geom_tile(data = join, aes(x = x, y = y, fill = `JSDM - SDM`)) +
    scale_fill_gradient2(na.value = "darkblue") +
    geom_contour(data = recla_sdm, aes(x = x, y = y, z = th_pred),
                 breaks = c(0.5), colour = "black", linewidth = 0.4) +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_bw() +
    coord_equal() +
    facet_grid(short_species ~ rcp) +
    theme(legend.position = "bottom")
  
  ggsave("./figures/diff_figure_V3.jpg", height = 6, width = 8)
}


combine.slope.map <- function(past_slope, single_past_slope) {
  
  
  label_data <- past_slope %>%
    group_by(variable) %>%
    summarise(
      x = mean(range(x, na.rm = TRUE)),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y = y_min + 0.17 * (y_max - y_min)
    )
  
  jsdm_plot <- ggplot() +
    geom_tile(data = filter(past_slope, is_gray),
              aes(x = x, y = y), fill = "grey80") +
    geom_tile(data = filter(past_slope, !is_gray),
              aes(x = x, y = y, fill = color_value)) +
    facet_wrap(~ variable) +
    scale_fill_gradient2(
      name = "Suit. trend\n(Slope)",
      low = "red", mid = "lightyellow", high = "darkgreen",
      midpoint = 0, breaks = c(-0.005, 0, 0.005),
      na.value = "darkblue"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.02))) +
    coord_equal() +
    theme_minimal() +
    xlab("Longitude") + ylab("Latitude") +
    theme(
      panel.spacing = unit(0.5, "pt"),   
      strip.text = element_blank(),     
      strip.background = element_blank()
    ) +
    geom_label(
      data = label_data,
      aes(x = x, y = y, label = variable),
      inherit.aes = FALSE,
      size = 3.5,
      label.padding = unit(2, "pt"),
      label.size = 0.25,
      fill = "white",
      alpha = 0.9,
      vjust = 1     
    )
  
  
  sdm_plot <- ggplot() +
    geom_tile(data = filter(single_past_slope, is_gray),
              aes(x = x, y = y), fill = "grey80") +
    geom_tile(data = filter(single_past_slope, !is_gray),
              aes(x = x, y = y, fill = color_value)) +
    facet_wrap(~ variable) +
    scale_fill_gradient2(
      name = "Suit. trend\n(Slope)",
      low = "red", mid = "lightyellow", high = "darkgreen",
      midpoint = 0, breaks = c(-0.005, 0, 0.005),
      na.value = "darkblue"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.02))) +
    coord_equal() +
    theme_minimal() +
    xlab("Longitude") + ylab("Latitude") +
    theme(
      panel.spacing = unit(0.5, "pt"),  
      strip.text = element_blank(),      
      strip.background = element_blank()
    ) +
    geom_label(
      data = label_data,
      aes(x = x, y = y, label = variable),
      inherit.aes = FALSE,
      size = 3.5,
      label.padding = unit(2, "pt"),
      label.size = 0.25,
      fill = "white",
      alpha = 0.9,
      vjust = 1     
    )
  
  arr <- ggarrange(
    sdm_plot,        
    jsdm_plot,       
    ncol = 1,
    heights = c(1, 1),   
    labels = c("a", "b"),
    common.legend = TRUE,
    legend = "bottom") 
  
  
  ggsave(
    filename = "./figures/combined_slope_v4.png",
    plot = arr,
    width = 7, height = 8, units = "in",
    dpi = 300,
    device = ragg::agg_png,
    bg = "white"
  )
}


plot.diff.sup <- function(fut_subset, single_fut_subset, th, myspecies) {

  jsdm <- fut_subset %>%
    mutate(jsdm = value) %>%
    select(-value)

  sdm <- single_fut_subset %>%
    mutate(sdm = value) %>%
    select(-value)

  join <- full_join(jsdm, sdm) %>%
    mutate("JSDM - SDM" = jsdm - sdm) %>%
    mutate(short_species = case_when(
      variable == "Abies alba" ~ "A. alba",
      variable == "Abies pinsapo" ~ "A. pinsapo",
      variable == "Cedrus atlantica" ~ "C. atlantica",
      variable == "Pinus halepensis" ~ "P. halepensis",
      variable == "Pinus nigra" ~ "P. nigra",
      variable == "Pinus pinea" ~ "P. pinea",
      variable == "Pinus sylvestris" ~ "P. sylvestris",
      variable == "Quercus coccifera" ~ "Q. coccifera",
      variable == "Quercus faginea" ~ "Q. faginea",
      variable == "Quercus ilex" ~ "Q. ilex",
      variable == "Quercus pyrenaica" ~ "Q. pyrenaica",
      variable == "Quercus rotundifolia" ~ "Q. rotundifolia")) %>%
    mutate(rcp = case_when(
      rcp == "rcp2.6" ~ "RCP 2.6",
      rcp == "rcp4.5" ~ "RCP 4.5",
      rcp == "rcp6.0" ~ "RCP 6.0",
      rcp == "rcp8.5" ~ "RCP 8.5",
      .default = rcp
    ))


  recla_sdm <- th.pred(th[[2]], single_fut_subset, myspecies) %>%
    mutate(rcp = case_when(
      rcp == "rcp2.6" ~ "RCP 2.6",
      rcp == "rcp4.5" ~ "RCP 4.5",
      rcp == "rcp6.0" ~ "RCP 6.0",
      rcp == "rcp8.5" ~ "RCP 8.5",
      .default = rcp
    ))



  ggplot() +
    geom_tile(data = join, aes(x = x, y = y, fill = `JSDM - SDM`)) +
    scale_fill_gradient2(na.value = "darkblue") +
    geom_contour(data = recla_sdm, aes(x = x, y = y, z = th_pred),
                 breaks = c(0.5), colour = "black", linewidth = 0.4) +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_bw() +
    coord_equal() +
    facet_grid(short_species ~ rcp) +
    theme(legend.position = "bottom")

  ggsave("./figures/diff_sup_v2.jpg", height = 14, width = 8)
}