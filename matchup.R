#matchup.R
#Author: Emanuele Spirito @ CNR-IREA
#Contributions: Gabriele Candiani, Lorenzo Parigi, Monica Pepe, Mirco Boschetti, Vittorio Brando
#Year: 2025
#Purpose: create matchup between hypstar station and satellite data

#get_hypstar_data ----
get_hypstar_data <- function(date_start,
                             date_end,
                             home_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/PROCESSED",
                             bands_HYP = get_bands_station("HYPSTAR")){
  library(tidytable)
  
  folder_list <- data.table::data.table(folder_name = list.dirs(home_path, recursive = FALSE, full.names = FALSE)) %>%
    mutate(date = substr(folder_name, 4, 11) %>% as.Date(., format = "%Y%m%d"),
           time = substr(folder_name, 12, 20)) %>%
    filter(date >= date_start,
           date <= date_end)
  
  res <- lapply(1:nrow(folder_list), function(row_n){
    folder_name <- folder_list$folder_name[row_n]
    date <- folder_list$date[row_n] 
    time <- folder_list$time[row_n] 
    
    file_list <- list.files(paste0(home_path, "/", folder_name), pattern = "L2A_REF", 
                            full.names = TRUE)
    
    if (length(file_list) == 0) {
      return(NULL)
    }
    
    ncin <- ncdf4::nc_open(file_list)
    
    hyp_ref <- ncdf4::ncvar_get(ncin, "reflectance")
    hyp_std <- ncdf4::ncvar_get(ncin, "std_reflectance")
    hyp_wvl  <-ncdf4::ncvar_get(ncin, "wavelength")
    hyp_vaa <-ncdf4::ncvar_get(ncin, "viewing_azimuth_angle")
    hyp_paa   <-ncdf4::ncvar_get(ncin, "pointing_azimuth_angle")
    hyp_vza  <-ncdf4::ncvar_get(ncin, "viewing_zenith_angle")
    hyp_sza <- ncdf4::ncvar_get(ncin, "solar_zenith_angle")
    hyp_saa <- ncdf4::ncvar_get(ncin, "solar_azimuth_angle")
    
    complete_ref <- data.table::data.table(date = date,
                                          time = time,
                                          VZA = hyp_vza,
                                          VAA = hyp_vaa,
                                          PAA = hyp_paa,
                                          SZA = hyp_sza,
                                          SAA = hyp_saa) %>%
      bind_cols(., hyp_ref %>% 
                  data.table::data.table() %>% 
                  `names<-`(paste0("B", 1:1554)) %>%
                  mutate(across(everything(), ~ as.integer(.x*10000)))) %>%
      mutate(across(c(VZA, VAA, PAA, SZA, SAA), ~round(.x)))
    
    complete_ref_long <- pivoting_longer_function(df = complete_ref, 
                             bands = bands_HYP) 
    
    complete_std <- data.table::data.table(date = date,
                                           time = time,
                                           VZA = hyp_vza,
                                           VAA = hyp_vaa,
                                           PAA = hyp_paa,
                                           SZA = hyp_sza,
                                           SAA = hyp_saa) %>%
      bind_cols(., hyp_std %>% 
                  data.table::data.table() %>% 
                  `names<-`(paste0("B", 1:1554)) %>%
                  mutate(across(everything(), ~ as.integer(.x*10000)))) %>%
      mutate(across(c(VZA, VAA, PAA, SZA, SAA), ~round(.x)))
    
    complete_std_long <- pivoting_longer_function(df = complete_std, 
                                                  bands = bands_HYP) %>%
      rename(std = Ref)
    
    complete_DB <- complete_ref_long %>% left_join(complete_std_long)
    
    return(complete_DB)
    
  }) %>% bind_rows() 
  
  # max value @ 836 nm < 6000
  filtered <- res %>%
    mutate(validate_B1000 = ifelse(bands == "B1000" & Ref > 6000,1,-1)) %>% 
    group_by(date,time,VZA,VAA,SAA,SZA,PAA) %>%
    mutate(max_validate_B1000 = max(validate_B1000)) %>%
    ungroup() %>%
    filter(max_validate_B1000 < 0) %>%
    select(-max_validate_B1000,-validate_B1000)
  
  #scale data
  scaled <- filtered %>% mutate(Ref = Ref/10000, std = std/10000)
  
  #remove PAA and convert time
  output <- scaled  %>% select(-PAA) %>%
    mutate(time = substring(time,2,8)) %>%
    mutate(time = sub("\\s+$", "", gsub('(.{2})', '\\1:', time))) %>%
    mutate(time = hms::as_hms(substring(time,1,8)))
  
  return(output)
}

#hypstar_spectral_resampling ----
hypstar_spectral_resampling <- function(dataset_res, #station data
                                        sat_sensor, #which satellite sensor to resample to
                                        bands_HYP = get_bands_station("HYPSTAR")){
  library(tidytable)
  
  bands_sensor <- get_bands_sat(sat_sensor)
  
  if(sat_sensor == "EnMAP"){
    sensor_mode <- "FWHM"
    which_bands <- paste0("B",1:166)
  }
  if(sat_sensor == "PRISMA"){
    sensor_mode <- "FWHM"
    which_bands <- paste0("B",1:151)
  }
  if(sat_sensor %in% c("S2A","S2B","S2C")){
    sensor_mode <- "SRF"
    which_bands <- c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12")
  }
  if(!(sat_sensor %in% c("PRISMA","EnMAP","S2A","S2B","S2C"))){
    print("ERRORE in spectral resampling")
  }
  
  source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/F01_Spectral_Resampling_V02.R")
  
  origin_band_centers <- sort(unique(bands_HYP$wvl))
  
  #Crea la matrice di ricampionamento
  SRF_sat <- suppressWarnings(SRF_adapter(origin_band_centers = origin_band_centers,
                                            mode = as.character(sensor_mode),
                                            sensor = as.character(sat_sensor),
                                            ref_folder = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/config_data")[,which_bands])
  
  #Eseguo il ricampionamento
  unpivoted_res <- pivoting_wider_function(unpivoted_with_bands = T,
                          output_pivoted = F,
                          df = dataset_res)
  
  band_variables <- names(unpivoted_res)[grepl(x = names(unpivoted_res), pattern = "^B.*")]
  
  file_resampled <- dataframe_resampler(total_df = unpivoted_res, 
                                        band_variables = band_variables, #qui serve il nome che le bande hanno nel file originale
                                        SRF = SRF_sat) 
  
  file_pivoted <- pivoting_longer_function(df = file_resampled,
                                             bands = bands_sensor)
  
  return(file_pivoted)
  
}

#get_s2_data ----
get_s2_data <- function(date_start, 
                        date_end,
                        hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                        repo_root_path = "//10.0.1.243/Repository/nr_products/sentinel2/l2/T32TQQ/" #masked S2 data
                        ){
  library(tidytable)
  year_start <- lubridate::year(as.Date(date_start))
  year_end <- lubridate::year(as.Date(date_end))
  years <- as.character(year_start:year_end)
  repo_year_path <- paste0(repo_root_path, years, "/")
  
  #now you enter in the folder where there are products
  get_products_paths <- function(year_path){
    list_of_products <- list.files(year_path, full.names = T)
    list_of_dates <- as.Date(substring(basename(list_of_products),8,15),"%Y%m%d")
    return(list_of_products[list_of_dates > date_start & list_of_dates < date_end])
  }
  
  product_paths <- unlist(lapply(repo_year_path, get_products_paths))
  
  hypstar_point <- sf::st_read(hypstar_point_path)
  
  get_four_closest_pixels_statistic <- function(band){
    #print(band)
    #create a raster with a buffer of 10 m around the hypstar point
    hypstar_point_buffer <- sf::st_buffer(hypstar_point,dist = 10)#10 is for 4 closest pixels
    #sf::st_write(hypstar_point_buffer,"//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/hypstar_buffer.shp")
    band_raster <- invisible(terra::rast(band))
    
    four_closest_pixels <- NULL
    four_closest_pixels <- tryCatch(t(terra::extract(band_raster, hypstar_point_buffer, fun = "table", 
                                          touches = F, ID = F, weights = T)), 
                                    error = function(e){NULL})
    if(is.null(four_closest_pixels)){
      output <- NA
    }else{
      value <- as.numeric(row.names(four_closest_pixels))
      
      four_closest_pixels <- data.table::as.data.table(four_closest_pixels)
      
      four_closest_pixels$value <- value
      
      names(four_closest_pixels) <- c("weight","value")
      
      four_closest_pixels <- four_closest_pixels[order(weight, decreasing = T)][1:4]#4 is 4 closest pixels
      
      #we do not return pixels because there are grid problems
      
      #statistics to return
      output <- c(mean(four_closest_pixels$value, na.rm = T),sd(four_closest_pixels$value, na.rm = T))
    }
    
    return(output)
  }
  
  get_statistics_dataset <- function(product_path){
    library(tidytable)
    
    sensor_ID <- substring(basename(product_path),24,26)
    
    list_of_tifs <- list.files(product_path, full.names = T, recursive = T,
                               pattern = "\\.tif$")
    
    bands <- list_of_tifs[!stringr::str_detect(list_of_tifs, "mask")]
    
    extract_names <- stringr::str_match(bands, "m_\\s*(.*?)\\s*.tif")[,2]
    
    bands_statistics <- lapply(bands, get_four_closest_pixels_statistic)
    
    bands_statistics_mean <- data.table::as.data.table(bands_statistics)[1,]
    
    bands_statistics_sd <- data.table::as.data.table(bands_statistics)[2,]
    
    statistics_dataset <- data.table::data.table()
    
    angles <- statistics_dataset %>% mutate(bands = extract_names) %>%
      mutate(value = as.numeric(bands_statistics_mean)) %>%
      filter(!grep(bands, pattern = "^B.*")) %>%
      pivot_wider(names_from = bands, values_from = value)
    
    statistics_pivoted <- statistics_dataset %>% mutate(bands = extract_names) %>%
      mutate(Ref = as.numeric(bands_statistics_mean)) %>%
      mutate(sd = as.numeric(bands_statistics_sd)) %>%
      filter(grep(bands, pattern = "^B.*")) %>%
      mutate(date = as.Date(substring(basename(product_path),8,15),"%Y%m%d")) %>%
      mutate(time = substring(basename(product_path),16,22)) %>% 
      mutate(time = substring(time,2,8)) %>%
      mutate(time = sub("\\s+$", "", gsub('(.{2})', '\\1:', time))) %>%
      mutate(time = hms::as_hms(substring(time,1,8))) %>%
      #pivot_wider(names_from = bands, values_from = Ref) %>%
      mutate(sensor = sensor_ID) %>%
      bind_cols(angles) %>%
      data.table::as.data.table()
    
    #statistics_dataset$bands <- extract_names
    
    #statistics_dataset$Ref <- as.numeric(bands_statistics)
    
    #here we want to extract also SZA, SAA, VZA, VAA
    
    return(statistics_pivoted)
  }
  
  all_datasets <- lapply(product_paths,get_statistics_dataset)
  
  binded <- dplyr::bind_rows(all_datasets)%>%
    na.omit() 
  
  with_wvl <- add_wvl(binded) %>%
    na.omit() 
  
  return(with_wvl)
}


#get_PRISMA_data ----
get_PRISMA_data <- function(date_start, 
                            date_end,
                            PRISMA_level_product = "L1",
                            processing_level = "Atcor_regrid_crop_smooth",
                            hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                            repo_root_path = "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/" #PRISMA data
                            ){
  #here we need to choose the format of PRISMA data
  library(tidytable)
  year_start <- lubridate::year(as.Date(date_start))
  year_end <- lubridate::year(as.Date(date_end))
  years <- as.character(year_start:year_end)
  repo_year_path <- paste0(repo_root_path, years, "/", PRISMA_level_product, "/")
  
  #now you enter in the folder where there are products
  get_products_paths <- function(year_path){
    list_of_products <- list.files(year_path, full.names = T)
    if(grepl("OFFL",basename(list_of_products[1]))){
      OFFL_offset <- 5
    }else{
      OFFL_offset <- 0
    }
    if(grepl("STD",basename(list_of_products[1]))){
      STD_offset <- 4
    }else{
      STD_offset <- 0
    }
    list_of_dates <- as.Date(substring(basename(list_of_products),6+OFFL_offset+STD_offset+nchar(PRISMA_level_product),13+OFFL_offset+STD_offset+nchar(PRISMA_level_product)),"%Y%m%d")
    return(list_of_products[list_of_dates > date_start & list_of_dates < date_end])
  }
  
  product_paths <- unlist(lapply(repo_year_path, get_products_paths))
  
  hypstar_point <- sf::st_read(hypstar_point_path)
  
  get_four_closest_pixels_statistic <- function(band){
    #print(band)
    #create a raster with a buffer of 10 m around the hypstar point
    hypstar_point_buffer <- sf::st_buffer(hypstar_point,dist = 30)#30 is for 4 closest pixels
    #sf::st_write(hypstar_point_buffer,"//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/hypstar_buffer.shp")
    #band_raster <- invisible(terra::rast(band))
    
    four_closest_pixels <- NULL
    four_closest_pixels <- tryCatch(t(terra::extract(band, hypstar_point_buffer, fun = "table", 
                                                     touches = F, ID = F, weights = T)), 
                                    error = function(e){NULL})
    if(is.null(four_closest_pixels)){
      output <- NA
    }else{
      value <- as.numeric(row.names(four_closest_pixels))
      
      four_closest_pixels <- data.table::as.data.table(four_closest_pixels)
      
      four_closest_pixels$value <- value
      
      names(four_closest_pixels) <- c("weight","value")
      
      four_closest_pixels <- four_closest_pixels[order(weight, decreasing = T)][1:4]#4 is 4 closest pixels
      
      #we do not return pixels because there are grid problems
      
      #statistics to return
      output <- c(mean(four_closest_pixels$value, na.rm = T),sd(four_closest_pixels$value, na.rm = T))
    }
    
    return(output)
  }
  
  get_statistics_dataset <- function(product_path){
    library(tidytable)
    
    print(product_path)
    
    sensor_ID <- "PRISMA"
    
    image_path <- list.files(paste0(product_path, "/", processing_level, "/"), full.names = T, recursive = T,
                               pattern = "\\.tif$")
    
    image <- terra::rast(image_path)
    
    create_list_of_rasters <- function(band){
      return(image[[band]])
    }
    
    bands <- lapply(1:230, create_list_of_rasters)
    
    bands_statistics <- suppressWarnings(lapply(bands, get_four_closest_pixels_statistic))
    
    bands_statistics_mean <- data.table::as.data.table(bands_statistics)[1,]
    
    bands_statistics_sd <- data.table::as.data.table(bands_statistics)[2,]
    
    statistics_dataset <- data.table::data.table()
    
    read_angle_and_time <- read.csv(paste0(product_path,"/ATCOR/all_angles_file.csv"))
    
    statistics_pivoted <- statistics_dataset %>% mutate(bands = paste0("B",1:230)) %>%
      mutate(Ref = as.numeric(bands_statistics_mean)) %>%
      mutate(sd = as.numeric(bands_statistics_sd)) %>%
      mutate(date = as.Date(substring(read_angle_and_time$date, 1,10),"%Y-%m-%d")) %>%
      mutate(time = hms::as_hms(substring(read_angle_and_time$date, 12,19))) %>% 
      #pivot_wider(names_from = bands, values_from = c(Ref,sd)) %>%
      mutate(sensor = sensor_ID) %>%
      mutate(SAA = read_angle_and_time$sunaz) %>%
      mutate(SZA = read_angle_and_time$sunzen) %>%
      mutate(VZA = read_angle_and_time$sensor_zenith) %>%
      mutate(VAA = read_angle_and_time$sensor_azimuth) %>%
      data.table::as.data.table()
    
    return(statistics_pivoted)
  }
  
  all_datasets <- lapply(product_paths,get_statistics_dataset)
  
  binded <- dplyr::bind_rows(all_datasets)%>%
    na.omit() %>% 
    mutate(Ref = Ref/10000) %>%
    mutate(sd = sd/10000)
  
  with_wvl <- add_wvl(binded) %>%
    na.omit() 
  
  return(with_wvl)
  
}

#get_EnMAP_data ----
get_EnMAP_data <- function(date_start, 
                            date_end,
                           EnMAP_level_product = "L2A",
                            hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                            repo_root_path = "//10.0.1.243/nr_data/3_rs_data/EnMAP/JDS/", #PRISMA data
                            bands_EnMAP = get_bands_sat("EnMAP")){
  #here we need to choose the format of PRISMA data
  library(tidytable)
  year_start <- lubridate::year(as.Date(date_start))
  year_end <- lubridate::year(as.Date(date_end))
  years <- as.character(year_start:year_end)
  repo_year_path <- paste0(repo_root_path, years, "/", EnMAP_level_product, "/")
  
  #now you enter in the folder where there are products
  get_products_paths <- function(year_path){
    list_of_products <- list.files(year_path, full.names = T)
    list_of_dates <- as.Date(substring(basename(list_of_products),30,37),"%Y%m%d")
    return(list_of_products[list_of_dates > date_start & list_of_dates < date_end])
  }
  
  product_paths <- unlist(lapply(repo_year_path, get_products_paths))
  
  hypstar_point <- sf::st_read(hypstar_point_path)
  
  get_four_closest_pixels_statistic <- function(band){
    #print(band)
    #create a raster with a buffer of 10 m around the hypstar point
    hypstar_point_buffer <- sf::st_buffer(hypstar_point,dist = 30)#30 is for 4 closest pixels
    #sf::st_write(hypstar_point_buffer,"//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/hypstar_buffer.shp")
    #band_raster <- invisible(terra::rast(band))
    
    four_closest_pixels <- NULL
    four_closest_pixels <- tryCatch(t(terra::extract(band, hypstar_point_buffer, fun = "table", 
                                                     touches = F, ID = F, weights = T)), 
                                    error = function(e){NULL})
    if(is.null(four_closest_pixels)){
      output <- NA
    }else{
      value <- as.numeric(row.names(four_closest_pixels))
      
      four_closest_pixels <- data.table::as.data.table(four_closest_pixels)
      
      four_closest_pixels$value <- value
      
      names(four_closest_pixels) <- c("weight","value")
      
      four_closest_pixels <- four_closest_pixels[order(weight, decreasing = T)][1:4]#4 is 4 closest pixels
      
      #we do not return pixels because there are grid problems
      
      #statistics to return
      output <- c(mean(four_closest_pixels$value, na.rm = T),sd(four_closest_pixels$value, na.rm = T))
    }
    
    return(output)
  }
  
  get_statistics_dataset <- function(product_path){
    library(tidytable)
    
    print(product_path)
    
    sensor_ID <- "EnMAP"
    
    image_path <- list.files(paste0(product_path, "/"), full.names = T, recursive = T,
                             pattern = "\\.BSQ$")
    
    metadata_path <- list.files(paste0(product_path, "/"), full.names = T, recursive = T,
                                pattern = "\\METADATA.XML$")
    
    #metadata_path <- "//10.0.1.243/nr_data/3_rs_data/EnMAP/JDS/2024/L2A/ENMAP01-____L2A-DT0000081830_20240714T104028Z_001_V010501_20241111T160454Z/ENMAP01-____L2A-DT0000081830_20240714T104028Z_001_V010501_20241111T160454Z-METADATA.XML"
    
    image <- terra::rast(image_path)
    
    create_list_of_rasters <- function(band){
      return(image[[band]])
    }
    
    bands <- lapply(1:224, create_list_of_rasters)
    
    bands_statistics <- suppressWarnings(lapply(bands, get_four_closest_pixels_statistic))
    
    bands_statistics_mean <- data.table::as.data.table(bands_statistics)[1,]
    
    bands_statistics_sd <- data.table::as.data.table(bands_statistics)[2,]
    
    data <- XML::xmlParse(metadata_path)
    
    xml_data <- XML::xmlToList(data)
    
    sunElevationAngle <- as.numeric(xml_data[["specific"]][["sunElevationAngle"]][["center"]][["text"]])
    if(sunElevationAngle < 0){
      print("Warning!. Elevation Angle less than zero. ")
    }
    sunZenithAngle <- 90-sunElevationAngle
    sunAzimithAngle <- as.numeric(xml_data[["specific"]][["sunAzimuthAngle"]][["center"]][["text"]])
    viewingZenithAngle <- as.numeric(xml_data[["specific"]][["viewingZenithAngle"]][["center"]][["text"]])
    viewingAzimithAngle <- as.numeric(xml_data[["specific"]][["viewingAzimuthAngle"]][["center"]][["text"]])
    
    statistics_dataset <- data.table::data.table()
    
    statistics_pivoted <- statistics_dataset %>% mutate(bands = paste0("B",1:224)) %>%
      mutate(Ref = as.numeric(bands_statistics_mean)) %>%
      mutate(sd = as.numeric(bands_statistics_sd)) %>%
      mutate(date = as.Date(substring(basename(product_path),30,37),"%Y%m%d")) %>%
      mutate(time = substring(basename(product_path),38,44)) %>% 
      mutate(time = substring(time,2,8)) %>%
      mutate(time = sub("\\s+$", "", gsub('(.{2})', '\\1:', time))) %>%
      mutate(time = hms::as_hms(substring(time,1,8))) %>%
      mutate(sensor = sensor_ID) %>%
      mutate(SAA = sunAzimithAngle) %>%
      mutate(SZA = sunZenithAngle) %>%
      mutate(VZA = viewingZenithAngle) %>%
      mutate(VAA = viewingAzimithAngle) %>%
      data.table::as.data.table()
    
    #statistics_dataset$bands <- extract_names
    
    #statistics_dataset$Ref <- as.numeric(bands_statistics)
    
    #here we want to extract also SZA, SAA, VZA, VAA
    
    return(statistics_pivoted)
  }
  
  all_datasets <- lapply(product_paths,get_statistics_dataset)
  
  binded <- dplyr::bind_rows(all_datasets)
  
  with_wvl <- add_wvl(binded) %>%
    na.omit() 
  
  return(with_wvl)
}

#matchup_function ----
matchup_function <- function(sat_data, 
                                      hypstar_data,
                                      delVZA = c(""), #vector deleting View Zenith Angle inside hypstar dataset
                                      delVAA = c(""), #vector deleting View Azimuth Angle inside hypstar datasets
                                      
                                      only_same_date = T,#only same date comparison accepted (true) or closest date (false)
                                      
                                      time_or_sun_angle = F, #consider time (true) or over sun azimuth and zenith angle (false) for filter
                                      just_closest_time_or_sun_angle = T,#closest time/angle (true) or all times/angles in a range (false)
                                      #weighted_mean = F,#to execute a mean (T) or a weighted mean (F) relative to the time/angle difference
                                      range_time_or_sun_angle = NULL,#if just_closest = F, here we need a range of hours/angles. -range/2,x,+range/2
                                      
                                      just_closest_view_angle = T,#closest time/angle (true) or all times/angles in a range (false)
                                      #weighted_mean = F,#to execute a mean (T) or a weighted mean (F) relative to the time/angle difference
                                      range_view_angle = NULL,#if just_closest = F, here we need a range of hours/angles. -range/2,x,+range/2
                                      
                             bands_HYP = get_bands_station("HYPSTAR")){
  
  #to remove specific sensor view angles
  hyp_pivot <- hypstar_data %>%     
    filter(!(VZA %in% delVZA)) %>%
    filter(!(VAA %in% delVAA))
  
  sat_name <- unique(sat_data$sensor)
  
  if(sat_name == "EnMAP"){
    hypstar_resampled <- hypstar_spectral_resampling(sat_sensor = "EnMAP",
                                                     dataset_res = hyp_pivot) %>% mutate(sensor = sat_name)
  }
  if(sat_name == "PRISMA"){
    hypstar_resampled <- hypstar_spectral_resampling(sat_sensor = "PRISMA",
                                                     dataset_res = hyp_pivot) %>% mutate(sensor = sat_name)
  }
  if(sat_name %in% c("S2A","S2B","S2C")){
    
    hypstar_resampled_to_S2A <- hypstar_spectral_resampling(sat_sensor = "S2A",
                                                     dataset_res = hyp_pivot) %>%
      mutate(sensor = "S2A")
    
    hypstar_resampled_to_S2B <- hypstar_spectral_resampling(sat_sensor = "S2B",
                                                            dataset_res = hyp_pivot)  %>%
      mutate(sensor = "S2B")
    
    hypstar_resampled_to_S2C <- hypstar_spectral_resampling(sat_sensor = "S2C",
                                                            dataset_res = hyp_pivot)  %>%
      mutate(sensor = "S2C")
    
    hypstar_resampled <- hypstar_resampled_to_S2A %>% bind_rows(hypstar_resampled_to_S2B) %>%
      bind_rows(hypstar_resampled_to_S2C)
  }
  if(!(sat_name %in% c("EnMAP","PRISMA","S2A","S2B","S2C"))){
    print("ERRORE")
  }
  
  bands_sat <- get_bands_sat(sat_name)
  
  #change col names not to confuse EnMAP B# and hypstar B#
  sat_col_renamed <- sat_data %>% rename_with(~paste0(.,"_sat"),-c(bands,date,sensor,wvl)) 
  
  hypstar_col_renamed <- hypstar_resampled %>% rename_with(~paste0(.,"_hyp"),-c(bands,date,sensor,wvl)) 
  
  
  library(tidytable)
  
  
  #DATE PROBLEM
  date_output <- date_function(sat_col_renamed = sat_col_renamed,
                               hypstar_col_renamed = hypstar_col_renamed,
                               only_same_date = only_same_date)
  
  #SUN TIME OR SUN ANGLE PROBLEM
  
  time_or_sun_angle_output <- time_or_sun_angle_function(date_output = date_output,
                                                         time_or_sun_angle = time_or_sun_angle,
                                                         just_closest_time_or_sun_angle = just_closest_time_or_sun_angle,
                                                         range_time_or_sun_angle = range_time_or_sun_angle)
  
  #VIEW ANGLE PROBLEM
  view_angle_output <- view_angle_function(time_or_sun_angle_output = time_or_sun_angle_output,
                                           range_view_angle = range_view_angle,
                                           just_closest_view_angle = just_closest_view_angle) 
  
  return(view_angle_output)
  
}


#view_mean ----
view_mean <- function(med_view_angle_output, range_view_angle){
  #filtro sugli angoli
  filtered <- med_view_angle_output %>%
    filter(angle_diff < range_view_angle/2) %>%
    select(-angle_diff) 
    
  ref_averaged <- filtered %>%
    group_by(date,bands,wvl,sensor) %>%
    select(-tos_sd_Ref_hyp) %>% #in futuro dobbiamo considerarla
    mutate(v_mean_Ref_hyp = mean(tos_mean_Ref_hyp)) %>%
    mutate(v_sd_Ref_hyp = sd(tos_mean_Ref_hyp)) %>%
    select(-tos_mean_Ref_hyp) #questa serve per fare la media e possiamo eliminarla
  
  valid_angle <- ref_averaged %>%
    mutate(count = row_number()) %>%
    mutate(count_max = max(count)) %>%
    select(-count) %>%
    mutate(validation_for_VZA_hyp = abs(max(VZA_hyp) - min(VZA_hyp)) > 0) %>%
    mutate(validation_for_VAA_hyp = abs(max(VAA_hyp) - min(VAA_hyp)) > 0) %>% 
    mutate(validation = validation_for_VZA_hyp & validation_for_VAA_hyp & count_max > 2) %>%
    select(-validation_for_VZA_hyp,-validation_for_VAA_hyp,-count_max)
  
  valid_angle_T <- valid_angle %>%
    filter(validation) %>%
    select(-VZA_hyp) %>%
    select(-VAA_hyp) %>%
    ungroup() %>% 
    unique()
  
  if(nrow(valid_angle_T) > 0){
    valid_angle_T <- valid_angle %>%
      filter(validation) %>%
      mutate(order = contoureR::orderPoints(cbind(VAA_hyp,VZA_hyp)[,1], cbind(VAA_hyp,VZA_hyp)[,2])) %>%
      arrange(order) %>% 
      select(-order) %>%
      mutate(VAA_hyp = ifelse(VAA_hyp > 180, VAA_hyp-360, VAA_hyp )) %>%
      mutate(mean_VZA_hyp = geosphere::centroid(cbind(VAA_hyp,VZA_hyp))[2]) %>%
      mutate(mean_VAA_hyp = geosphere::centroid(cbind(VAA_hyp,VZA_hyp))[1]) %>%
      mutate(mean_VAA_hyp = ifelse(mean_VAA_hyp < 0, 360 + mean_VAA_hyp, mean_VAA_hyp)) %>%
      ungroup() %>%
      select(-VZA_hyp) %>%
      select(-VAA_hyp) %>%
      unique() 
  }
  
  valid_angle_F <- valid_angle %>%
    filter(!validation) %>%
    select(-VZA_hyp) %>%
    select(-VAA_hyp) %>%
    ungroup() %>% 
    unique()
  
  if(nrow(valid_angle_F) > 0){
    valid_angle_F <- valid_angle %>%
      filter(!validation) %>%
      mutate(mean_VZA_hyp = mean(VZA_hyp)) %>%
      mutate(mean_VAA_hyp = mean(VAA_hyp)) %>%
      ungroup() %>%
      select(-VZA_hyp) %>%
      select(-VAA_hyp)%>%
      unique() 
  }
  
  
  view_angle_output <- valid_angle_T %>% bind_rows(valid_angle_F) %>%
    select(-validation)
  
  return(view_angle_output)
}

#view_closest ----
view_closest <- function(med_output){
  cols_sd <- c(sd_hyp = "NA")
  
  cols <- c(operations = "NA")
  
  write_ops <- "closest_view_angle"
  
  closest_output <- med_output %>% 
      group_by(date,bands,wvl,sensor) %>%
      mutate(lowest_angle_diff = min(angle_diff)) %>%
      ungroup() %>%
      filter(angle_diff == lowest_angle_diff) %>%
      select(-angle_diff,-lowest_angle_diff) 
  
  closest_output_check <- closest_output %>%
    tibble::add_column(!!!cols_sd[!names(cols_sd) %in% names(.)]) %>%
    mutate(sd_hyp = ifelse(is.na(sd_hyp),0,sd_hyp)) %>% #when the time is the closest there is no standard deviation
    tibble::add_column(!!!cols[!names(cols) %in% names(.)]) %>%
    mutate(operations = paste(operations,write_ops,sep=";") ) 
  
  return(closest_output_check)
}

#time_or_sun_closest ----
time_or_sun_closest <- function(med_output, time_or_angle_diff){
  
  cols_sd <- c(sd_hyp = "NA")
    
  cols <- c(operations = "NA")
  
  if(time_or_angle_diff){
    write_ops <- "closest_time"
    closest_output <- med_output %>% 
      group_by(date,bands,sensor,wvl) %>%
      mutate(lowest_time_diff = min(time_diff)) %>%
      ungroup() %>%
      filter(time_diff == lowest_time_diff) %>%
      select(-time_diff,-lowest_time_diff) 
  }else{
    write_ops <- "closest_sun_angle"
    closest_output <- med_output %>% 
      group_by(date,bands,wvl,sensor) %>%
      mutate(lowest_angle_diff = min(angle_diff)) %>%
      ungroup() %>%
      filter(angle_diff == lowest_angle_diff) %>%
      select(-angle_diff,-lowest_angle_diff) 
  }
  
  closest_output_check <- closest_output %>%
    tibble::add_column(!!!cols_sd[!names(cols_sd) %in% names(.)]) %>%
    mutate(sd_hyp = ifelse(is.na(sd_hyp),0,sd_hyp)) %>% #when the time is the closest there is no standard deviation
    group_by(date,sensor,bands,wvl) %>% #qui c'è il problema che scegliendo l'orario più vicino ci possono essere due sun angles differenti che ci portiamo dietro 
    mutate(SZA_hyp = mean(SZA_hyp)) %>% 
    mutate(SAA_hyp = mean(SAA_hyp)) %>%
    ungroup() %>%
    unique() %>% 
    tibble::add_column(!!!cols[!names(cols) %in% names(.)]) %>%
    mutate(operations = paste(operations,write_ops,sep=";") ) 
  
  return(closest_output_check)
}

#time_or_sun_mean ----
time_or_sun_mean <- function(med_output, range_time_or_sun_angle, time_or_sun_angle){
  
  cols_sd <- c(sd_hyp = NA_real_)
  
  cols <- c(operations = "")
  #filtro
  if(time_or_sun_angle){
    write_ops <- paste0("average_in_time_window_of_",range_time_or_sun_angle,"hours")
    filtered <- med_output %>%
      filter(time_diff < range_time_or_sun_angle/2*3600) %>%
      select(-time_diff) 
  }else{
    write_ops <- paste0("average_in_sun_angle_window_of",range_time_or_sun_angle,"degrees")
    filtered <- med_output %>%
      filter(angle_diff < range_time_or_sun_angle/2) %>%
      select(-angle_diff) 
  }
  
  ref_averaged <- filtered %>%
    group_by(date,bands,sensor,wvl,VZA_hyp,VAA_hyp) %>%
    mutate(count = row_number()) %>%
    mutate(count_max = max(count)) %>%
    select(-count) %>%
    mutate(mean_Ref_hyp = mean(Ref_hyp)) %>%
    tibble::add_column(!!!cols_sd[!names(cols_sd) %in% names(.)]) %>%
    #mutate(sd_hyp = case_match(sd_hyp,NA ~ 0)) %>% #poi rimetto questo
    mutate(sd_hyp = jitter(0,1))%>%
    mutate(mean_sd_hyp = error_propation(sd_hyp)) %>%
    select(-Ref_hyp,-sd_hyp) %>%
    rename(Ref_hyp = mean_Ref_hyp, sd_hyp = mean_sd_hyp) %>%
    mutate(time_hyp = mean(time_hyp))
  
  valid_angle <- ref_averaged %>%
    mutate(validation_for_SZA_hyp = abs(max(SZA_hyp) - min(SZA_hyp)) > 0) %>%
    mutate(validation_for_SAA_hyp = abs(max(SAA_hyp) - min(SAA_hyp)) > 0) %>% 
    mutate(validation = validation_for_SZA_hyp & validation_for_SAA_hyp & count_max > 2) %>%
    select(-validation_for_SZA_hyp,-validation_for_SAA_hyp,-count_max) 
  
  valid_angle_T <- valid_angle %>%
    filter(validation) %>% 
    ungroup() %>% 
    unique()
  
  if(nrow(valid_angle_T)>0){
    valid_angle_T <- valid_angle %>%
      filter(validation) %>%
      mutate(order = contoureR::orderPoints(cbind(SAA_hyp,SZA_hyp)[,1], cbind(SAA_hyp,SZA_hyp)[,2])) %>%
      arrange(order) %>% 
      select(-order) %>%
      mutate(SAA_hyp = ifelse(SAA_hyp > 180, SAA_hyp-360, SAA_hyp )) %>%
      mutate(mean_SZA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[2]) %>%
      mutate(mean_SAA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[1]) %>%
      mutate(mean_SAA_hyp = ifelse(mean_SAA_hyp < 0, 360 + mean_SAA_hyp, mean_SAA_hyp)) %>%
      ungroup() %>%
      select(-SZA_hyp) %>%
      select(-SAA_hyp) %>%
      rename(SZA_hyp = mean_SZA_hyp) %>%
      rename(SAA_hyp = mean_SAA_hyp) %>%
      unique() 
  }
  
  valid_angle_F <- valid_angle %>%
    filter(!validation) %>% 
    ungroup() %>% 
    unique()
  
  if(nrow(valid_angle_F)>0){
    valid_angle_F <- valid_angle %>%
      filter(!validation) %>%
      mutate(mean_SZA_hyp = mean(SZA_hyp)) %>%
      mutate(mean_SAA_hyp = mean(SAA_hyp)) %>%
      ungroup() %>%
      select(-SZA_hyp) %>%
      select(-SAA_hyp) %>%
      rename(SZA_hyp = mean_SZA_hyp) %>%
      rename(SAA_hyp = mean_SAA_hyp) %>%
      unique() 
  }
  
  time_or_sun_angle_output <- valid_angle_T %>% bind_rows(valid_angle_F) %>%
    select(-validation) %>%
    tibble::add_column(!!!cols[!names(cols) %in% names(.)]) %>%
    mutate(operations = paste(operations,write_ops,sep=";") ) %>% 
    group_by(date,sensor,wvl,bands) %>%#there could be a little difference in the mean since some images could not pass quality checks so we make another average over angles
    mutate(SAA_hyp = mean(SAA_hyp)) %>%
    mutate(SZA_hyp = mean(SZA_hyp)) %>%
    ungroup() %>%
    unique()
  
  return(time_or_sun_angle_output)
}
#date_function ----
date_function <- function(sat_col_renamed, hypstar_col_renamed, only_same_date){
  if(only_same_date){
    #only_same_date
    date_output <- sat_col_renamed %>% left_join(hypstar_col_renamed, by = c("bands","date","sensor","wvl")) %>%
      na.omit()
    
  }else{
    #closest date 
    get_name <- unique(sat_col_renamed$sensor)
    print(paste0("closest date will have the ",get_name," date as 'date' variable"))
    sat_with_date <- sat_col_renamed %>% rename(date_sat = date)
    hypstar_with_date <- hypstar_col_renamed %>% rename(date_hyp = date)
    date_output <- sat_with_date %>% left_join(hypstar_with_date, by = c("bands","sensor","wvl")) %>%
      mutate(date_diff = abs(difftime(date_sat, date_hyp))) %>%
      group_by(sensor,bands) %>%
      mutate(lowest_date_diff = min(date_diff)) %>%
      ungroup() %>%
      filter(date_diff == lowest_date_diff) %>%
      select(-date_diff,-lowest_date_diff) %>%
      select(-date_hyp) %>%
      rename(date = date_sat) %>%
      na.omit()
  }
  
  return(date_output)
}

#time_or_sun_angle_function ----
time_or_sun_angle_function <- function(date_output,time_or_sun_angle,just_closest_time_or_sun_angle, range_time_or_sun_angle){
  cols <- c(operations = "")
  
  if(time_or_sun_angle){
    
    #create time variables and time difference
    time_output <- date_output %>%
      mutate(time_diff = abs(difftime(time_sat, time_hyp, units = "sec")))
    
    if(just_closest_time_or_sun_angle){
      #only closest time
      time_or_sun_angle_output <- time_or_sun_closest(time_output, 
                                                      time_or_angle_diff = T) 
      
    }else{
      #take a window of times and mean
      
      time_or_sun_angle_output <- time_or_sun_mean(time_output,range_time_or_sun_angle,time_or_sun_angle) 
      
      
      
    }
  }else{
    #create time variables and angle difference
    sun_angle_output <- date_output %>%
      mutate(angle_diff = 180/pi*acos( sin(SAA_sat*pi/180)*sin(SAA_hyp*pi/180) + cos(SAA_sat*pi/180) * cos(SAA_hyp*pi/180) * cos( SZA_sat*pi/180-SZA_hyp*pi/180 ))) 
    #https://en.wikipedia.org/wiki/Great-circle_distance#cite_ref-2
    
    if(just_closest_time_or_sun_angle){
      #only closest angle
      time_or_sun_angle_output <- time_or_sun_closest(sun_angle_output, 
                                                      time_or_angle_diff = F) %>%
        rename(tos_mean_Ref_hyp = Ref_hyp) %>%#QUIIIIIII controllato
        mutate(tos_sd_Ref_hyp = 0)%>%
        rename(mean_time_hyp = time_hyp)
    }else{
      #take a window of angles and mean
      
      time_or_sun_angle_output <- time_or_sun_mean(sun_angle_output,range_time_or_sun_angle,time_or_sun_angle)
      
    }
  }
  return(time_or_sun_angle_output)
}

#view_angle_function ----
view_angle_function <- function(time_or_sun_angle_output, just_closest_view_angle, range_view_angle){
  med_view_angle_output <- time_or_sun_angle_output %>%
    mutate(angle_diff = 180/pi*acos( sin(VAA_sat*pi/180)*sin(VAA_hyp*pi/180) + cos(VAA_sat*pi/180) * cos(VAA_hyp*pi/180) * cos( VZA_sat*pi/180-VZA_hyp*pi/180 ))) 
  #https://en.wikipedia.org/wiki/Great-circle_distance#cite_ref-2
  
  if(just_closest_view_angle){
    #only closest angle
    view_angle_output <- view_closest(med_view_angle_output, time_or_angle_diff = F)
    
  }else{
    #take a window of angles and mean
    
    view_angle_output <- view_mean(med_view_angle_output, range_view_angle)
    
  }
  return(view_angle_output)
}

#pivoting_longer_function ----
add_wvl <- function(df){
  library(tidytable)
  
  sat_name <- unique(df$sensor)
  
  sat_bands <- get_bands_sat(sat_name)
  
  res <- tidytable::left_join(df,sat_bands) %>% filter(!is.na(Ref))
  
  return(res)
}

#pivoting_longer_function ----
pivoting_longer_function <- function(df, bands){
  res2_long <- tidyr::pivot_longer(df, cols = contains("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
  res2_long <- tidytable::left_join(res2_long,bands) %>% filter(!is.na(Ref)) #not to use na.omit() that could be dangerous
  
  return(res2_long)
}

#pivoting_wider_function ----
pivoting_wider_function <- function(output_pivoted, df, unpivoted_with_bands){
  library(tidytable)
  
  if(!output_pivoted){
    if(unpivoted_with_bands){
      not_pivoted <- df %>% select(-wvl) %>%
        tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
    }else{
      not_pivoted <- df %>% select(-bands) %>%
        tidyr::pivot_wider(names_from = "wvl", values_from = "Ref")
    }
    
  }else{
    not_pivoted <- df
  }
  
  return(not_pivoted)
}

#get_bands_sat ----
get_bands_sat <- function(sat_name){
  if(sat_name == "PRISMA"){
    bands_csv <- read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/bands_wavelenght/PRISMA_bands.csv")
  }
  if(sat_name == "EnMAP"){
    bands_csv <- read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/bands_wavelenght/EnMAP_bands.csv")
  }
  if(sat_name %in% c("S2A","S2B","S2C","S2")){
    bands_csv <- read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/bands_wavelenght/s2_bands.csv")
  }
  if(!(sat_name %in% c("S2A","S2B","S2C","S2","PRISMA","EnMAP"))){
    bands_csv <- NULL
  }
  return(bands_csv)
}

#get_bands_station ----
get_bands_station <- function(station_name){
  if(station_name == "HYPSTAR"){
    bands_csv <- read.csv("//10.0.1.243/projects/2023_HyperNet-POP/6-script/bands_HYP.csv")
  }else{
    bands_csv <- NULL
  }
  return(bands_csv)
}

#plot_ribbon ----
plot_ribbon <- function(matchup_result,
                         color_station = "darkgreen",
                         color_sat_EnMAP = "blue",
                         color_sat_PRISMA = "orange",
                         color_sat_S2 = "red",
                         alpha = .1,
                         title = "Titolo di prova",
                         save_path = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/last_plot.png"){
  library(ggplot2)
  
  sat_name <- unique(matchup_result$sensor)
  
  if(sat_name == "PRISMA"){
    sat_color <- color_sat_PRISMA
  }
  if(sat_name == "EnMAP"){
    sat_color <- color_sat_EnMAP
  }
  if(sat_name %in% c("S2A","S2B","S2C","S2")){
    sat_color <- color_sat_S2
  }
  
  plot <- ggplot(data = matchup_result) + 
    geom_line(aes(x = wvl, y = Ref_sat), linewidth = 1, color = sat_color) +
    geom_line(aes(x = wvl, y = v_mean_Ref_hyp), linewidth = 1, color = color_station) +
    geom_ribbon(aes(x = wvl, y = Ref_sat, ymin = Ref_sat - sd_sat, ymax = Ref_sat + sd_sat), alpha = alpha, fill = sat_color) +
    geom_ribbon(aes(x = wvl, y = v_mean_Ref_hyp, ymin = v_mean_Ref_hyp - v_sd_Ref_hyp, ymax = v_mean_Ref_hyp + v_sd_Ref_hyp), alpha = alpha, fill = color_station) +
    ggtitle(title) +
    ylab("R [-]") +
    xlab("Wavelength [nm]")
  
  ggsave(plot = plot, filename = save_path, height = 3.5, width = 7)
  
  return(plot)
}

plot_rainbow <- function(matchup_result,
                         point_size = 3,
                         title = "Titolo di prova",
                         save_path = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/last_plot.png"){
  library(ggplot2)
  palette <- c(rev(rainbow(8)[-length(rainbow(8))]),"#bbbbbb","#777777","#444444","#000000")
  values_palette <- scales::rescale(x = c(350,400,450,500,550,600,650,700,1000,1300,1800),
                            to = c(0,1))
  sat_name <- unique(matchup_result$sensor)
  
  
  
  plot <- ggplot() +
    geom_point(data = matchup_result, 
               mapping = aes(x = v_mean_Ref_hyp, 
                             y = Ref_sat,
                             color = wvl), 
               size = point_size)+
    scale_color_gradientn(colours = palette,
                          na.value = NA,
                          values = values_palette,
                          name = "wavelenghts [nm]")+
    ggtitle(title) +
    ylab(paste0(sat_name," R [-]")) +
    geom_abline(slope = 1,
                intercept = 0,
                alpha = 0.4,
                linetype = "dashed")+
    xlab("Hypstar R [-]")
  
  ggsave(plot = plot, filename = save_path, height = 3.5, width = 7)
  
  return(plot)
}

plot_relative_errors <- function(matchup_result,
                                 point_size = 3,
                                 title = "Titolo di prova",
                                 save_path = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/last_plot.png"){
  relative_error_df <- matchup_result %>% mutate(relative_error = (Ref_sat-v_mean_Ref_hyp)/ v_mean_Ref_hyp *100) 
  library(ggplot2)
  plot <- ggplot(data = relative_error_df) +
    geom_abline(slope = 0, 
                intercept = 0,
                linetype = "dashed") +
    geom_point(mapping = aes(x = wvl, y = relative_error, color = relative_error))+
    scale_color_gradient2('Relative error (%)',      # titolo della legenda
                          low = 'red', mid = 'darkgreen', high = 'red',
                          midpoint = 0)+
    ylab("Relative error (%)")+
    xlab("wavelengths [nm]")
    
  ggsave(plot = plot, filename = save_path, height = 3.5, width = 7)
  
  return(plot)
}

error_propation <- function(sd_hyp){
  length <- length(sd_hyp)
  squares <- as.numeric(lapply(sd_hyp,function(x){return(x^2)}))
  sum <- sum(squares)
  output <- sqrt(sum)/length
  return(output)
}
