#matchup.R
#Author: Emanuele Spirito @ CNR-IREA
#Contributions: Gabriele Candiani, Lorenzo Parigi, Monica Pepe, Mirco Boschetti, Vittorio Brando
#Year: 2025
#Purpose: create matchup between hypstar station and satellite data

#get_hypstar_data ----
get_hypstar_data <- function(date_start,
                             date_end,
                             home_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/PROCESSED",
                             bands_HYP = get_bands_station("HYPSTAR_JDS"),
                             filename = NULL,
                             output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/"){
  #library(tidytable)
  
  folder_list <- data.table::data.table(folder_name = list.dirs(home_path, recursive = FALSE, full.names = FALSE)) |>
    tidytable::mutate(date = as.Date(substr(folder_name, 4, 11), format = "%Y%m%d"),
           time = substr(folder_name, 12, 20)) |>
    tidytable::filter(date >= date_start,
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
    hyp_u_rel_random <- ncdf4::ncvar_get(ncin, "u_rel_random_reflectance")
    hyp_u_rel_systematic <- ncdf4::ncvar_get(ncin, "u_rel_systematic_reflectance")
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
                                          SAA = hyp_saa) |>
      tidytable::bind_cols(hyp_ref |> 
                  data.table::data.table() |> 
                  `names<-`(paste0("B", 1:1554)) ) |>
      tidytable::mutate(tidytable::across(c(VZA, VAA, PAA, SZA, SAA), ~round(.x)))
    
    complete_ref_long <- pivoting_longer_function(df = complete_ref, 
                             bands = bands_HYP) 
    
    complete_std <- data.table::data.table(date = date,
                                           time = time,
                                           VZA = hyp_vza,
                                           VAA = hyp_vaa,
                                           PAA = hyp_paa,
                                           SZA = hyp_sza,
                                           SAA = hyp_saa) |>
      tidytable::bind_cols(hyp_std |> 
                  data.table::data.table() |> 
                  `names<-`(paste0("B", 1:1554)) ) |>
      tidytable::mutate(tidytable::across(c(VZA, VAA, PAA, SZA, SAA), ~round(.x)))
    
    complete_std_long <- pivoting_longer_function(df = complete_std, 
                                                  bands = bands_HYP) |>
      tidytable::rename(std = Ref)
    
    # complete_u_rel_random <- data.table::data.table(date = date,
    #                                        time = time,
    #                                        VZA = hyp_vza,
    #                                        VAA = hyp_vaa,
    #                                        PAA = hyp_paa,
    #                                        SZA = hyp_sza,
    #                                        SAA = hyp_saa) |>
    #   bind_cols(., hyp_u_rel_random |> 
    #               data.table::data.table() |> 
    #               `names<-`(paste0("B", 1:1554)) ) |>
    #   mutate(across(c(VZA, VAA, PAA, SZA, SAA), ~round(.x)))
    # 
    # complete_u_rel_random_long <- pivoting_longer_function(df = complete_u_rel_random, 
    #                                               bands = bands_HYP) |>
    #   rename(u_rel_random = Ref)
    # 
    # complete_u_rel_systematic <- data.table::data.table(date = date,
    #                                                 time = time,
    #                                                 VZA = hyp_vza,
    #                                                 VAA = hyp_vaa,
    #                                                 PAA = hyp_paa,
    #                                                 SZA = hyp_sza,
    #                                                 SAA = hyp_saa) |>
    #   bind_cols(., hyp_u_rel_systematic |> 
    #               data.table::data.table() |> 
    #               `names<-`(paste0("B", 1:1554)) ) |>
    #   mutate(across(c(VZA, VAA, PAA, SZA, SAA), ~round(.x)))
    # 
    # complete_u_rel_systematic_long <- pivoting_longer_function(df = complete_u_rel_systematic, 
    #                                                   bands = bands_HYP) |>
    #   rename(u_rel_systematic = Ref)
    # 
    complete_DB <- complete_ref_long |> tidytable::left_join(complete_std_long) #|>
      # left_join(complete_u_rel_random_long) |> left_join(complete_u_rel_systematic_long)
    
    return(complete_DB)
    
  }) |> tidytable::bind_rows() 
  
  # max value @ 836 nm < 6000
  filtered <- res |>
    tidytable::mutate(validate_B1000 = ifelse(bands == "B1000" & Ref > 6000,1,-1)) |> 
    tidytable::group_by(date,time,VZA,VAA,SAA,SZA,PAA) |>
    tidytable::mutate(max_validate_B1000 = max(validate_B1000)) |>
    tidytable::ungroup() |>
    tidytable::filter(max_validate_B1000 < 0) |>
    tidytable::select(-max_validate_B1000,-validate_B1000)
  
  #scale data
  # scaled <- filtered |> mutate(Ref = Ref/10000, std = std/10000,
  #                               u_rel_random = u_rel_random/10000,
  #                               u_rel_systematic = u_rel_systematic/10000)
  # 
  #remove PAA and convert time
  output <- filtered  |> tidytable::select(-PAA) |>
    tidytable::mutate(time = substring(time,2,8)) |>
    tidytable::mutate(time = sub("\\s+$", "", gsub('(.{2})', '\\1:', time))) |>
    tidytable::mutate(time = hms::as_hms(substring(time,1,8)))
  
  if(!is.null(filename)){
    save_function(output, filename, output_directory)
  }
  
  return(output)
}

#hypstar_spectral_resampling ----
hypstar_spectral_resampling <- function(dataset_res, #station data
                                        sat_sensor, #which satellite sensor to resample to
                                        bands_HYP = get_bands_station("HYPSTAR_JDS")){
  # library(tidytable)
  
  bands_sensor <- get_bands_sat(substr(sat_sensor,1,2))
  
  if(substr(sat_sensor,1,2) == "En"){
    sensor_mode <- "FWHM"
    which_bands <- paste0("B",1:166)
  }
  if(substr(sat_sensor,1,2) == "PR"){
    sensor_mode <- "FWHM"
    which_bands <- paste0("B",1:151)
  }
  if(substr(sat_sensor,1,2) == "S2"){
    sensor_mode <- "SRF"
    which_bands <- c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12")
  }
  if(!(substr(sat_sensor,1,2) %in% c("PR","En","S2"))){
    print("ERRORE in spectral resampling")
  }
  
  source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/F01_Spectral_Resampling_V02.R")
  # library(tidytable)
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
  
  # detach("package:tidytable", unload = TRUE)
  
  return(file_pivoted)
  
}

#get_s2_data ----
get_s2_data <- function(date_start, 
                        date_end,
                        hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                        repo_root_path = "//10.0.1.243/Repository/nr_products/sentinel2/l2/T32TQQ/", #masked S2 data
                        filename = NULL,
                        output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/",
                        scale = F){
  # library(tidytable)
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
    # library(tidytable)
    
    sensor_ID <- substring(basename(product_path),24,26)
    
    list_of_tifs <- list.files(product_path, full.names = T, recursive = T,
                               pattern = "\\.tif$")
    
    bands <- list_of_tifs[!stringr::str_detect(list_of_tifs, "mask")]
    
    extract_names <- stringr::str_match(bands, "m_\\s*(.*?)\\s*.tif")[,2]
    
    bands_statistics <- lapply(bands, get_four_closest_pixels_statistic)
    
    bands_statistics_mean <- data.table::as.data.table(bands_statistics)[1,]
    
    bands_statistics_sd <- data.table::as.data.table(bands_statistics)[2,]
    
    statistics_dataset <- data.table::data.table()
    
    angles <- statistics_dataset |> tidytable::mutate(bands = extract_names) |>
      tidytable::mutate(value = as.numeric(bands_statistics_mean)) |>
      tidytable::filter(!grep(bands, pattern = "^B.*")) |>
      tidytable::pivot_wider(names_from = bands, values_from = value)
    
    statistics_pivoted <- statistics_dataset |> tidytable::mutate(bands = extract_names) |>
      tidytable::mutate(Ref = as.numeric(bands_statistics_mean)) |>
      tidytable::mutate(sd = as.numeric(bands_statistics_sd)) |>
      tidytable::filter(grep(bands, pattern = "^B.*")) |>
      tidytable::mutate(date = as.Date(substring(basename(product_path),8,15),"%Y%m%d")) |>
      tidytable::mutate(time = substring(basename(product_path),16,22)) |> 
      tidytable::mutate(time = substring(time,2,8)) |>
      tidytable::mutate(time = sub("\\s+$", "", gsub('(.{2})', '\\1:', time))) |>
      tidytable::mutate(time = hms::as_hms(substring(time,1,8))) |>
      #pivot_wider(names_from = bands, values_from = Ref) |>
      tidytable::mutate(sensor = sensor_ID) |>
      tidytable::bind_cols(angles) |>
      data.table::as.data.table()
    
    #statistics_dataset$bands <- extract_names
    
    #statistics_dataset$Ref <- as.numeric(bands_statistics)
    
    #here we want to extract also SZA, SAA, VZA, VAA
    
    return(statistics_pivoted)
  }
  
  all_datasets <- lapply(product_paths,get_statistics_dataset)
  
  binded <- dplyr::bind_rows(all_datasets)|>
    stats::na.omit() 
  
  if(scale){
    binded <- binded |> 
      tidytable::mutate(Ref = Ref/10000) |>
      tidytable::mutate(sd = sd/10000)
  }
  
  with_wvl <- add_wvl(binded) |>
    stats::na.omit() 
  
  if(!is.null(filename)){
    save_function(with_wvl, filename, output_directory)
  }
  
  return(with_wvl)
}


#get_PRISMA_data ----
get_PRISMA_data <- function(date_start, 
                            date_end,
                            PRISMA_level_product = "L1",
                            processing_level = "Atcor_regrid_crop_smooth",
                            hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                            repo_root_path = "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/", #PRISMA data
                            filename = NULL,
                            output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/",
                            scale = F){
  #here we need to choose the format of PRISMA data
  # library(tidytable)
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
    # library(tidytable)
    
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
    
    statistics_pivoted <- statistics_dataset |> tidytable::mutate(bands = paste0("B",1:230)) |>
      tidytable::mutate(Ref = as.numeric(bands_statistics_mean)) |>
      tidytable::mutate(sd = as.numeric(bands_statistics_sd)) |>
      tidytable::mutate(date = as.Date(substring(read_angle_and_time$date, 1,10),"%Y-%m-%d")) |>
      tidytable::mutate(time = hms::as_hms(substring(read_angle_and_time$date, 12,19))) |> 
      #pivot_wider(names_from = bands, values_from = c(Ref,sd)) |>
      tidytable::mutate(sensor = sensor_ID) |>
      tidytable::mutate(SAA = read_angle_and_time$sunaz) |>
      tidytable::mutate(SZA = read_angle_and_time$sunzen) |>
      tidytable::mutate(VZA = read_angle_and_time$sensor_zenith) |>
      tidytable::mutate(VAA = read_angle_and_time$sensor_azimuth) |>
      data.table::as.data.table()
    
    return(statistics_pivoted)
  }
  
  all_datasets <- lapply(product_paths,get_statistics_dataset)
  
  binded <- dplyr::bind_rows(all_datasets)|>
    stats::na.omit()
  
  if(scale){
     binded <- binded |> 
       tidytable::mutate(Ref = Ref/10000) |>
       tidytable::mutate(sd = sd/10000)
  }
  
  with_wvl <- add_wvl(binded) |>
    stats::na.omit() 
  
  if(!is.null(filename)){
    save_function(with_wvl, filename, output_directory)
  }
  
  return(with_wvl)
  
}

#get_EnMAP_data ----
get_EnMAP_data <- function(date_start, 
                            date_end,
                           EnMAP_level_product = "L2A", #enmap level product
                            hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp", #hypstar shapefile path
                            repo_root_path = "//10.0.1.243/nr_data/3_rs_data/EnMAP/JDS/", #the main folder of PRISMA data
                            # bands_EnMAP = get_bands_sat(substr("EnMAP",1,2)), #where to find enmap bands
                           filename = NULL, #if not NULL the enmap data will be saved with this filename
                           output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/", #output directory if filename is not null
                           scale = F){
  # library(tidytable)
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
    # library(tidytable)
    
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
    
    statistics_pivoted <- statistics_dataset |> tidytable::mutate(bands = paste0("B",1:224)) |>
      tidytable::mutate(Ref = as.numeric(bands_statistics_mean)) |>
      tidytable::mutate(sd = as.numeric(bands_statistics_sd)) |>
      tidytable::mutate(date = as.Date(substring(basename(product_path),30,37),"%Y%m%d")) |>
      tidytable::mutate(time = substring(basename(product_path),38,44)) |> 
      tidytable::mutate(time = substring(time,2,8)) |>
      tidytable::mutate(time = sub("\\s+$", "", gsub('(.{2})', '\\1:', time))) |>
      tidytable::mutate(time = hms::as_hms(substring(time,1,8))) |>
      tidytable::mutate(sensor = sensor_ID) |>
      tidytable::mutate(SAA = sunAzimithAngle) |>
      tidytable::mutate(SZA = sunZenithAngle) |>
      tidytable::mutate(VZA = viewingZenithAngle) |>
      tidytable::mutate(VAA = viewingAzimithAngle) |>
      data.table::as.data.table()
    
    #statistics_dataset$bands <- extract_names
    
    #statistics_dataset$Ref <- as.numeric(bands_statistics)
    
    #here we want to extract also SZA, SAA, VZA, VAA
    
    return(statistics_pivoted)
  }
  
  all_datasets <- lapply(product_paths,get_statistics_dataset)
  
  binded <- dplyr::bind_rows(all_datasets)
  
  if(scale){
    binded <- binded |> 
      tidytable::mutate(Ref = Ref/10000) |>
      tidytable::mutate(sd = sd/10000)
  }
  
  with_wvl <- add_wvl(binded) |>
    stats::na.omit() 
  
  if(!is.null(filename)){
    save_function(with_wvl, filename, output_directory)
  }
  
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
                             range_time_or_sun_angle = NULL,#if just_closest = F, here we need a range of hours/angles. -range/2,x,+range/2
                             
                             just_closest_view_angle = T,#closest time/angle (true) or all times/angles in a range (false)
                             #weighted_mean = F,#to execute a mean (T) or a weighted mean (F) relative to the time/angle difference
                             range_view_angle = NULL,#if just_closest = F, here we need a range of hours/angles. -range/2,x,+range/2
                             filename = NULL, #a custom name for file to save, otherwise it will be given a standard file name
                             output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/", #where to save output of matchup (and polar plot)
                             # bands_HYP = get_bands_station("HYPSTAR_JDS"), #bands of Hypstar station
                             print_polar = F, # if we want to save in the output_directory a polar plot (T) or not (F)
                             weighted_mean = F #if we want a weighted mean (T) or a simple mean (F)
){
  
  #to remove specific sensor view angles
  hyp_pivot <- hypstar_data |>     
    tidytable::filter(!(VZA %in% delVZA)) |>
    tidytable::filter(!(VAA %in% delVAA))
  
  sat_sensor <- unique(substr(sat_data$sensor,1,2))
  
  if(sat_sensor == "En"){
    hypstar_resampled <- hypstar_spectral_resampling_ref_and_std(sat_sensor = "EnMAP",
                                                        dataset_res = hyp_pivot)
  }
  if(sat_sensor == "PR"){
    hypstar_resampled <- hypstar_spectral_resampling_ref_and_std(sat_sensor = "PRISMA",
                                                                 dataset_res = hyp_pivot)
  }
  if(sat_sensor == "S2"){
    hypstar_resampled_to_S2A <- hypstar_spectral_resampling_ref_and_std(sat_sensor = "S2A",
                                                                 dataset_res = hyp_pivot)

    hypstar_resampled_to_S2B <- hypstar_spectral_resampling_ref_and_std(sat_sensor = "S2B",
                                                                        dataset_res = hyp_pivot)
    
    hypstar_resampled_to_S2C <- hypstar_spectral_resampling_ref_and_std(sat_sensor = "S2C",
                                                                        dataset_res = hyp_pivot)
    
    hypstar_resampled <- hypstar_resampled_to_S2A |> tidytable::bind_rows(hypstar_resampled_to_S2B) |>
      tidytable::bind_rows(hypstar_resampled_to_S2C)
  }
  if(!(sat_sensor %in% c("En","PR","S2"))){
    print("ERRORE")
  }
  
  bands_sat <- get_bands_sat(sat_sensor)
  
  #change col names not to confuse EnMAP B# and hypstar B#
  sat_col_renamed <- sat_data |> tidytable::rename_with(~paste0(.,"_sat"),-c(bands,date,sensor,wvl)) 
  
  hypstar_col_renamed <- hypstar_resampled |> tidytable::rename_with(~paste0(.,"_hyp"),-c(bands,date,sensor,wvl)) 
  
  
  # library(tidytable)
  
  
  #DATE PROBLEM
  date_output <- date_function(sat_col_renamed = sat_col_renamed,
                               hypstar_col_renamed = hypstar_col_renamed,
                               only_same_date = only_same_date)
  
  #SUN TIME OR SUN ANGLE PROBLEM
  
  time_or_sun_angle_output <- time_or_sun_angle_function(date_output = date_output,
                                                         time_or_sun_angle = time_or_sun_angle,
                                                         just_closest_time_or_sun_angle = just_closest_time_or_sun_angle,
                                                         range_time_or_sun_angle = range_time_or_sun_angle,
                                                         weighted_mean = weighted_mean)
  
  if(print_polar){
    plot_polar(time_or_sun_angle_output = time_or_sun_angle_output,
               output_directory = gsub("csv","plots", output_directory))
  }
  
  #VIEW ANGLE PROBLEM
  view_angle_output <- view_angle_function(time_or_sun_angle_output = time_or_sun_angle_output,
                                           range_view_angle = range_view_angle,
                                           just_closest_view_angle = just_closest_view_angle,
                                           weighted_mean = weighted_mean) 
  
  #save
  if(!is.null(filename)){
    save_function(view_angle_output, filename, output_directory)
  }else{
    save_function(view_angle_output, paste0(gsub(";","_",unique(view_angle_output$operations)),".csv"), output_directory) 
  }
  
  
  return(view_angle_output)
  
}


#view_mean ----
view_mean <- function(med_view_angle_output, range_view_angle, weighted_mean){
  cols <- c(operations = "")
  cols_num <- c(number_data_in_window = "")
  
  write_ops <- paste0("average_in_view_angle_window_of",range_view_angle,"degrees")
  
  if(weighted_mean){
    write_ops <- paste0("weighted_",write_ops)
  }
  
  #filter over angles (keep only angles in the selected range)
  filtered <- med_view_angle_output |>
    tidytable::filter(angle_diff < range_view_angle/2) |>
    tidytable::mutate(weight = 1/abs(angle_diff)) |>
    tidytable::select(-angle_diff) 
    
  ref_averaged <- filtered |>
    tidytable::mutate(todoweight = weighted_mean) |>
    tidytable::mutate(weight = ifelse(todoweight,weight,1)) |>
    tidytable::group_by(date,bands,wvl,sensor) |>
    tidytable::mutate(mean_Ref_hyp = weighted.mean(x = Ref_hyp, w = weight)) |>
    tidytable::select(-weight,-todoweight) |>
    #mutate(mean_std_hyp = error_propagation(std_hyp)) |>
    #mutate(mean_std_hyp = sd(Ref_hyp)) |>
    tidytable::mutate(mean_std_hyp = error_propagation_v2(std_hyp = std_hyp,
                                               sample_std_hyp = sd(Ref_hyp)))
  
  names_cols_num <- names(ref_averaged)
  
  ref_averaged <- ref_averaged |>
    tibble::add_column(!!!cols_num[!names(cols_num) %in% names_cols_num]) |>
    tidytable::mutate(number_data_in_window = paste(number_data_in_window,length(std_hyp),sep=";")) |>
    tidytable::select(-Ref_hyp,-std_hyp) |>
    tidytable::rename(Ref_hyp = mean_Ref_hyp, 
           std_hyp = mean_std_hyp)
  
  if(!weighted_mean){
    #here we do not weight the mean because we want to find a value to express the center of the sky visualized from the data 
    valid_angle <- ref_averaged |>
      tidytable::mutate(count = tidytable::row_number()) |>
      tidytable::mutate(count_max = max(count)) |>
      tidytable::select(-count) |>
      tidytable::mutate(validation_for_VZA_hyp = abs(max(VZA_hyp) - min(VZA_hyp)) > 0) |>
      tidytable::mutate(validation_for_VAA_hyp = abs(max(VAA_hyp) - min(VAA_hyp)) > 0) |> 
      tidytable::mutate(validation = validation_for_VZA_hyp & validation_for_VAA_hyp & count_max > 2) |>
      tidytable::select(-validation_for_VZA_hyp,-validation_for_VAA_hyp,-count_max)
    
    #here we have to split the computation:
    #first case: the points in the sky are 3 or more but not aligned so that 
    #they form a polygon. For those we compute the centroids
    #second case: the points in the sky are 2 or more but aligned so that they form 
    #a segment. For those we compute the midpoint. 
    
    valid_angle_T <- valid_angle |>
      tidytable::filter(validation) |>
      tidytable::ungroup() |> 
      unique()
    
    if(nrow(valid_angle_T) > 0){
      valid_angle_T <- valid_angle |>
        tidytable::filter(validation) |>
        tidytable::mutate(order = inverse_order_points(VAA_hyp,VZA_hyp)) |> 
        tidytable::arrange(order) |> 
        tidytable::select(-order) |>
        tidytable::mutate(VAA_hyp = ifelse(VAA_hyp > 180, VAA_hyp-360, VAA_hyp )) |>
        tidytable::mutate(mean_VZA_hyp = geosphere::centroid(cbind(VAA_hyp,VZA_hyp))[2]) |>
        tidytable::mutate(mean_VAA_hyp = geosphere::centroid(cbind(VAA_hyp,VZA_hyp))[1]) |>
        tidytable::mutate(mean_VAA_hyp = ifelse(mean_VAA_hyp < 0, 360 + mean_VAA_hyp, mean_VAA_hyp)) |>
        tidytable::ungroup() |>
        tidytable::select(-VZA_hyp) |>
        tidytable::select(-VAA_hyp) |>
        tidytable::rename(VZA_hyp = mean_VZA_hyp) |>
        tidytable::rename(VAA_hyp = mean_VAA_hyp) |>
        unique() 
    }
    
    valid_angle_F <- valid_angle |>
      tidytable::filter(!validation) |>
      tidytable::ungroup() |> 
      unique()
    
    if(nrow(valid_angle_F) > 0){
      valid_angle_F <- valid_angle |>
        tidytable::filter(!validation) |>
        tidytable::mutate(VAA_hyp = ifelse(VAA_hyp > 180, VAA_hyp-360, VAA_hyp )) |>
        tidytable::arrange(VAA_hyp,VZA_hyp) |>
        tidytable::mutate(mean_VZA_hyp = geosphere::midPoint(cbind(VAA_hyp,VZA_hyp)[1,],cbind(VAA_hyp,VZA_hyp)[length(VAA_hyp),])[2]) |>
        tidytable::mutate(mean_VAA_hyp = geosphere::midPoint(cbind(VAA_hyp,VZA_hyp)[1,],cbind(VAA_hyp,VZA_hyp)[length(VAA_hyp),])[1]) |>
        tidytable::mutate(mean_VAA_hyp = ifelse(mean_VAA_hyp < 0, 360 + mean_VAA_hyp, mean_VAA_hyp)) |>
        tidytable::ungroup() |>
        tidytable::select(-VZA_hyp) |>
        tidytable::select(-VAA_hyp)|>
        tidytable::rename(VZA_hyp = mean_VZA_hyp) |>
        tidytable::rename(VAA_hyp = mean_VAA_hyp) |>
        unique() 
    }
    
    view_angle_output <- valid_angle_T |> tidytable::bind_rows(valid_angle_F) |>
      tidytable::select(-validation)
    
    names_cols <- names(view_angle_output)
    
    view_angle_output <- view_angle_output |>
      tibble::add_column(!!!cols[!names(cols) %in% names_cols]) |>
      tidytable::mutate(operations = paste(operations,write_ops,sep=";") )
    
  }else{
    #as requested from GC when we do the weighted mean the viewing angle geometry recorded is the one
    #of the satellite
    
    names_cols <- names(ref_averaged)
    
    view_angle_output <- ref_averaged |>
      tibble::add_column(!!!cols[!names(cols) %in% names_cols]) |>
      tidytable::mutate(operations = paste(operations,write_ops,sep=";") ) |> 
      tidytable::mutate(VAA_hyp = VAA_sat) |>
      tidytable::mutate(VZA_hyp = VZA_sat) |>
      tidytable::ungroup() |>
      unique()
  }
  
  
  return(view_angle_output)
}

#view_closest ----
view_closest <- function(med_view_angle_output){
  
  cols <- c(operations = "")
  
  write_ops <- "closest_view_angle"
  
  closest_output <- med_view_angle_output |> 
    tidytable::group_by(date,bands,wvl,sensor) |>
    tidytable::mutate(lowest_angle_diff = min(angle_diff)) |>
    tidytable::ungroup() |>
    tidytable::filter(angle_diff == lowest_angle_diff) |>
    tidytable::select(-angle_diff,-lowest_angle_diff) 
  
  names_cols <- names(closest_output)
  
  closest_output_check <- closest_output |>
    tibble::add_column(!!!cols[!names(cols) %in% names_cols]) |>
    tidytable::mutate(operations = paste(operations,write_ops,sep=";") ) 
  
  return(closest_output_check)
}

#time_or_sun_closest ----
time_or_sun_closest <- function(med_output, time_or_angle_diff){
  
  cols <- c(operations = "")
  
  if(time_or_angle_diff){
    write_ops <- "closest_time"
    closest_output <- med_output |> 
      tidytable::group_by(date,bands,sensor,wvl) |>
      tidytable::mutate(lowest_time_diff = min(time_diff)) |>
      tidytable::ungroup() |>
      tidytable::filter(time_diff == lowest_time_diff) |>
      tidytable::select(-time_diff,-lowest_time_diff) 
  }else{
    write_ops <- "closest_sun_angle"
    closest_output <- med_output |> 
      tidytable::group_by(date,bands,wvl,sensor) |>
      tidytable::mutate(lowest_angle_diff = min(angle_diff)) |>
      tidytable::ungroup() |>
      tidytable::filter(angle_diff == lowest_angle_diff) |>
      tidytable::select(-angle_diff,-lowest_angle_diff) 
  }
  
  closest_output_check <- closest_output |>
    #qui c'è il problema che scegliendo l'orario più vicino ci possono essere due sun angles differenti che ci portiamo dietro 
    tidytable::group_by(date,sensor,bands,wvl) |> 
    tidytable::mutate(SZA_hyp = mean(SZA_hyp)) |> 
    tidytable::mutate(SAA_hyp = mean(SAA_hyp)) |>
    tidytable::ungroup() |>
    unique() 
  
  names_cols <- names(closest_output_check)
  
  closest_output_check <- closest_output_check |>
    tibble::add_column(!!!cols[!names(cols) %in% names_cols]) |>
    tidytable::mutate(operations = paste(operations,write_ops,sep=";") ) 
  
  return(closest_output_check)
}

#time_or_sun_mean ----
time_or_sun_mean <- function(med_output, range_time_or_sun_angle, time_or_sun_angle, weighted_mean){
  
  cols <- c(operations = "")
  cols_num <- c(number_data_in_window = "")
  
  #filter over time/angle (keep only times/angles in the range)
  if(time_or_sun_angle){
    write_ops <- paste0("average_in_time_window_of_",range_time_or_sun_angle,"hours")
    filtered <- med_output |>
      tidytable::filter(time_diff < range_time_or_sun_angle/2*3600) |>
      tidytable::mutate(weight = 1/as.numeric(abs(time_diff))) |>
      tidytable::select(-time_diff) 
  }else{
    write_ops <- paste0("average_in_sun_angle_window_of",range_time_or_sun_angle,"degrees")
    filtered <- med_output |>
      tidytable::filter(angle_diff < range_time_or_sun_angle/2) |>
      tidytable::mutate(weight = 1/as.numeric(abs(angle_diff))) |>
      tidytable::select(-angle_diff)
  }
  
  if(weighted_mean){
    write_ops <- paste0("weighted_",write_ops)
  }
  
  ref_averaged <- filtered |>
    tidytable::mutate(todoweight = weighted_mean) |>
    tidytable::mutate(weight = ifelse(todoweight,weight,1)) |>
    tidytable::group_by(date,bands,sensor,wvl,VZA_hyp,VAA_hyp) |>
    tidytable::mutate(mean_Ref_hyp = weighted.mean(x = Ref_hyp, w = weight)) |>
    tidytable::select(-weight,-todoweight) |>
    #mutate(mean_std_hyp = error_propagation(std_hyp)) |>
    #mutate(mean_std_hyp = sd(Ref_hyp)) |>
    tidytable::mutate(mean_std_hyp = error_propagation_v2(std_hyp = std_hyp,
                                               sample_std_hyp = sd(Ref_hyp)))
  
  names_cols_num <- names(ref_averaged)
  
  ref_averaged <- ref_averaged |>
    tibble::add_column(!!!cols_num[!names(cols_num) %in% names_cols_num]) |>
    tidytable::mutate(number_data_in_window = paste(number_data_in_window,length(std_hyp),sep=";")) |>
    tidytable::select(-Ref_hyp,-std_hyp) |>
    tidytable::rename(Ref_hyp = mean_Ref_hyp, 
           std_hyp = mean_std_hyp) |>
    tidytable::mutate(time_hyp = mean(time_hyp))
  
  if(!weighted_mean){
    #here we do not weight the mean because we want to find a value to express the center of the sky visualized from the data 
    valid_angle <- ref_averaged |>
      tidytable::mutate(count = tidytable::row_number()) |>
      tidytable::mutate(count_max = max(count)) |>
      tidytable::select(-count) |>
      tidytable::mutate(validation_for_SZA_hyp = abs(max(SZA_hyp) - min(SZA_hyp)) > 0) |>
      tidytable::mutate(validation_for_SAA_hyp = abs(max(SAA_hyp) - min(SAA_hyp)) > 0) |> 
      tidytable::mutate(validation = validation_for_SZA_hyp & validation_for_SAA_hyp & count_max > 2) |>
      tidytable::select(-validation_for_SZA_hyp,-validation_for_SAA_hyp,-count_max) 
    
    #here we have to split the computation:
    #first case: the points in the sky are 3 or more but not aligned so that 
    #they form a polygon. For those we compute the centroids
    #second case: the points in the sky are 2 or more but aligned so that they form 
    #a segment. For those we compute the midpoint. 
    
    valid_angle_T <- valid_angle |>
      tidytable::filter(validation) |> 
      tidytable::ungroup() |> 
      unique()
    
    if(nrow(valid_angle_T)>0){
      valid_angle_T <- valid_angle |>
        tidytable::filter(validation) |>
        tidytable::mutate(order = inverse_order_points(SAA_hyp,SZA_hyp)) |>
        tidytable::arrange(order) |> 
        tidytable::select(-order) |>
        tidytable::mutate(SAA_hyp = ifelse(SAA_hyp > 180, SAA_hyp-360, SAA_hyp )) |>
        tidytable::mutate(mean_SZA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[2]) |>
        tidytable::mutate(mean_SAA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[1]) |>
        tidytable::mutate(mean_SAA_hyp = ifelse(mean_SAA_hyp < 0, 360 + mean_SAA_hyp, mean_SAA_hyp)) |>
        tidytable::ungroup() |>
        tidytable::select(-SZA_hyp) |>
        tidytable::select(-SAA_hyp) |>
        tidytable::rename(SZA_hyp = mean_SZA_hyp) |>
        tidytable::rename(SAA_hyp = mean_SAA_hyp) |>
        unique() 
    }
    
    valid_angle_F <- valid_angle |>
      tidytable::filter(!validation) |> 
      tidytable::ungroup() |> 
      unique()
    
    if(nrow(valid_angle_F)>0){
      valid_angle_F <- valid_angle |>
        tidytable::filter(!validation) |>
        tidytable::mutate(SAA_hyp = ifelse(SAA_hyp > 180, SAA_hyp-360, SAA_hyp )) |>
        tidytable::arrange(SAA_hyp,SZA_hyp) |>
        tidytable::mutate(mean_SZA_hyp = geosphere::midPoint(cbind(SAA_hyp,SZA_hyp)[1,],cbind(SAA_hyp,SZA_hyp)[length(SAA_hyp),])[2]) |>
        tidytable::mutate(mean_SAA_hyp = geosphere::midPoint(cbind(SAA_hyp,SZA_hyp)[1,],cbind(SAA_hyp,SZA_hyp)[length(SAA_hyp),])[1]) |>
        tidytable::mutate(mean_SAA_hyp = ifelse(mean_SAA_hyp < 0, 360 + mean_SAA_hyp, mean_SAA_hyp)) |>
        tidytable::ungroup() |>
        tidytable::select(-SZA_hyp) |>
        tidytable::select(-SAA_hyp) |>
        tidytable::rename(SZA_hyp = mean_SZA_hyp) |>
        tidytable::rename(SAA_hyp = mean_SAA_hyp) |>
        unique() 
    }
    
    time_or_sun_angle_output <- valid_angle_T |> tidytable::bind_rows(valid_angle_F) |>
      tidytable::select(-validation)
      
    names_cols <- names(time_or_sun_angle_output)
    
    #while different VZA and VAA are taken the sun position changes so its SAA and 
    #SZA change. This means that creating averages grouping by VZA and VAA 
    #will have different centroids in SAA and SZA. For this reason there will
    #be slightly different SZA and SAA for each combination of VZA and VAA.
    #To take a uniform value for all of them we average over all of them
    
    time_or_sun_angle_output <- time_or_sun_angle_output |>
      tibble::add_column(!!!cols[!names(cols) %in% names_cols]) |>
      tidytable::mutate(operations = paste(operations,write_ops,sep=";") ) |> 
      tidytable::group_by(date,sensor,wvl,bands) |>
      tidytable::mutate(SAA_hyp = mean(SAA_hyp)) |>
      tidytable::mutate(SZA_hyp = mean(SZA_hyp)) |>
      tidytable::ungroup() |>
      unique()
    
  }else{
    #as requested from GC when we do the weighted mean the viewing angle geometry recorded is the one
    #of the satellite
    
    names_cols <- names(ref_averaged)
    
    time_or_sun_angle_output <- ref_averaged |>
      tibble::add_column(!!!cols[!names(cols) %in% names_cols]) |>
      tidytable::mutate(operations = paste(operations,write_ops,sep=";") ) |> 
      tidytable::mutate(SZA_hyp = SZA_sat) |>
      tidytable::mutate(SAA_hyp = SAA_sat) |>
      tidytable::ungroup() |>
      unique()
  }
  
  return(time_or_sun_angle_output)
}

#date_function ----
date_function <- function(sat_col_renamed, hypstar_col_renamed, only_same_date){
  if(only_same_date){
    #only_same_date
    date_output <- sat_col_renamed |> tidytable::left_join(hypstar_col_renamed, by = c("bands","date","sensor","wvl")) |>
      stats::na.omit()
    
  }else{
    #closest date 
    get_name <- unique(sat_col_renamed$sensor)
    print(paste0("closest date will have the ",get_name," date as 'date' variable"))
    sat_with_date <- sat_col_renamed |> tidytable::rename(date_sat = date)
    hypstar_with_date <- hypstar_col_renamed |> tidytable::rename(date_hyp = date)
    date_output <- sat_with_date |> tidytable::left_join(hypstar_with_date, by = c("bands","sensor","wvl")) |>
      tidytable::mutate(date_diff = abs(difftime(date_sat, date_hyp))) |>
      tidytable::group_by(sensor,bands) |>
      tidytable::mutate(lowest_date_diff = min(date_diff)) |>
      tidytable::ungroup() |>
      tidytable::filter(date_diff == lowest_date_diff) |>
      tidytable::select(-date_diff,-lowest_date_diff) |>
      tidytable::select(-date_hyp) |>
      tidytable::rename(date = date_sat) |>
      stats::na.omit()
  }
  
  return(date_output)
}

#time_or_sun_angle_function ----
time_or_sun_angle_function <- function(date_output,time_or_sun_angle,just_closest_time_or_sun_angle, range_time_or_sun_angle, weighted_mean){
  cols <- c(operations = "")
  
  if(time_or_sun_angle){
    
    #create time variables and time difference
    time_output <- date_output |>
      tidytable::mutate(time_diff = abs(difftime(time_sat, time_hyp, units = "sec")))
    
    if(just_closest_time_or_sun_angle){
      #only closest time
      time_or_sun_angle_output <- time_or_sun_closest(med_output = time_output, 
                                                      time_or_angle_diff = T) 
    }else{
      #take a window of times and mean
      time_or_sun_angle_output <- time_or_sun_mean(med_output = time_output,
                                                   range_time_or_sun_angle = range_time_or_sun_angle,
                                                   time_or_sun_angle = time_or_sun_angle,
                                                   weighted_mean = weighted_mean) 
    }
  }else{
    #create time variables and angle difference
    sun_angle_output <- date_output |>
      tidytable::mutate(angle_diff = 180/pi*acos( sin(SAA_sat*pi/180)*sin(SAA_hyp*pi/180) + cos(SAA_sat*pi/180) * cos(SAA_hyp*pi/180) * cos( SZA_sat*pi/180-SZA_hyp*pi/180 ))) 
    #https://en.wikipedia.org/wiki/Great-circle_distance#cite_ref-2
    
    if(just_closest_time_or_sun_angle){
      #only closest angle
      time_or_sun_angle_output <- time_or_sun_closest(med_output = sun_angle_output, 
                                                      time_or_angle_diff = F) 
    }else{
      #take a window of angles and mean
      time_or_sun_angle_output <- time_or_sun_mean(med_output = sun_angle_output,
                                                   range_time_or_sun_angle = range_time_or_sun_angle,
                                                   time_or_sun_angle = time_or_sun_angle,
                                                   weighted_mean = weighted_mean) 
    }
  }
  return(time_or_sun_angle_output)
}

#view_angle_function ----
view_angle_function <- function(time_or_sun_angle_output, just_closest_view_angle, range_view_angle, weighted_mean){
  med_view_angle_output <- time_or_sun_angle_output |>
    tidytable::mutate(angle_diff = 180/pi*acos( sin(VAA_sat*pi/180)*sin(VAA_hyp*pi/180) + cos(VAA_sat*pi/180) * cos(VAA_hyp*pi/180) * cos( VZA_sat*pi/180-VZA_hyp*pi/180 ))) 
  #https://en.wikipedia.org/wiki/Great-circle_distance#cite_ref-2
  
  if(just_closest_view_angle){
    #only closest angle
    view_angle_output <- view_closest(med_view_angle_output = med_view_angle_output)
    
  }else{
    #take a window of angles and mean
    
    view_angle_output <- view_mean(med_view_angle_output = med_view_angle_output, 
                                   range_view_angle = range_view_angle,
                                   weighted_mean = weighted_mean)
    
  }
  return(view_angle_output)
}

#add_wvl ----
add_wvl <- function(df){
  # library(tidytable)
  
  sat_name <- unique(substr(df$sensor,1,2))
  
  sat_bands <- get_bands_sat(sat_name)
  
  res <- tidytable::left_join(df,sat_bands) |> tidytable::filter(!is.na(Ref))
  
  return(res)
}

#pivoting_longer_function ----
pivoting_longer_function <- function(df, bands){
  res2_long <- tidyr::pivot_longer(df, cols = contains("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
  res2_long <- tidytable::left_join(res2_long,bands) |> tidytable::filter(!is.na(Ref)) #not to use na.omit() that could be dangerous
  
  return(res2_long)
}

#pivoting_wider_function ----
pivoting_wider_function <- function(output_pivoted, df, unpivoted_with_bands){
  # library(tidytable)
  
  if(!output_pivoted){
    if(unpivoted_with_bands){
      not_pivoted <- df |> tidytable::select(-wvl) |>
        tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
    }else{
      not_pivoted <- df |> tidytable::select(-bands) |>
        tidyr::pivot_wider(names_from = "wvl", values_from = "Ref")
    }
    
  }else{
    not_pivoted <- df
  }
  
  return(not_pivoted)
}

#get_bands_sat ----
get_bands_sat <- function(sat_name){
  if(sat_name == "PR"){
    bands_csv <- read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/bands_wavelenght/PRISMA_bands.csv")
  }
  if(sat_name == "En"){
    bands_csv <- read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/bands_wavelenght/EnMAP_bands.csv")
  }
  if(sat_name == "S2"){
    bands_csv <- read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/bands_wavelenght/s2_bands.csv")
  }
  if(!(sat_name %in% c("S2","PR","En"))){
    bands_csv <- NULL
  }
  return(bands_csv)
}

#get_bands_station ----
get_bands_station <- function(station_name){
  if(station_name == "HYPSTAR_JDS"){
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
                         alpha = .3,
                        title = paste0(unique(substr(matchup_result$sensor,1,2)),"_vs_hypstar",unique(matchup_result$operations)),
                        output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/",
                        filename = paste0("ribbon_",gsub(";","_",title))){
  library(ggplot2)
  
  sat_name <- unique(substr(matchup_result$sensor,1,2))
  
  if(sat_name == "PR"){
    sat_color <- color_sat_PRISMA
  }
  if(sat_name == "En"){
    sat_color <- color_sat_EnMAP
  }
  if(sat_name == "S2"){
    sat_color <- color_sat_S2
  }
  
  plot <- ggplot(data = matchup_result) + 
    geom_line(aes(x = wvl, y = Ref_sat), linewidth = 1, color = sat_color) +
    geom_line(aes(x = wvl, y = Ref_hyp), linewidth = 1, color = color_station) +
    geom_ribbon(aes(x = wvl, y = Ref_sat, ymin = Ref_sat - sd_sat, ymax = Ref_sat + sd_sat), alpha = alpha, fill = sat_color) +
    geom_ribbon(aes(x = wvl, y = Ref_hyp, ymin = Ref_hyp - std_hyp, ymax = Ref_hyp + std_hyp), alpha = alpha, fill = color_station) +
    ggtitle(title) +
    ylab("R [-]") +
    xlab("Wavelength [nm]")
  
  save_path <- paste0(output_directory, filename,".png")
    
  ggsave(plot = plot, filename = save_path, height = 3.5, width = 7)
  
  return(plot)
}

#plot_rainbow ----
plot_rainbow <- function(matchup_result,
                         point_size = 3,
                         title = paste0(unique(substr(matchup_result$sensor,1,2)),"_vs_hypstar",unique(matchup_result$operations)),
                         output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/",
                         filename = paste0("rainbow_",gsub(";","_",title))){
  library(ggplot2)
  palette <- c(rev(rainbow(8)[-length(rainbow(8))]),"#bbbbbb","#777777","#444444","#000000")
  values_palette <- scales::rescale(x = c(350,400,450,500,550,600,650,700,1000,1300,1800),
                            to = c(0,1))
  sat_name <- unique(substr(matchup_result$sensor,1,2))
  
  
  
  plot <- ggplot() +
    geom_point(data = matchup_result, 
               mapping = aes(x = Ref_hyp, 
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
  
  save_path <- paste0(output_directory, filename,".png")
  
  ggsave(plot = plot, filename = save_path, height = 3.5, width = 7)
  
  return(plot)
}

#plot_relative_errors ----
plot_relative_errors <- function(matchup_result,
                                 point_size = 3,
                                 title = paste0(unique(substr(matchup_result$sensor,1,2)),"_vs_hypstar",unique(matchup_result$operations)),
                                 output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/",
                                 filename = paste0("relative_error_",gsub(";","_",title))){
  relative_error_df <- matchup_result |> tidytable::mutate(relative_error = (Ref_sat-Ref_hyp)/ Ref_hyp *100) 
  library(ggplot2)
  plot <- ggplot(data = relative_error_df) +
    geom_abline(slope = 0, 
                intercept = 0,
                linetype = "dashed") +
    ggtitle(title)+
    geom_point(mapping = aes(x = wvl, y = relative_error, color = relative_error))+
    scale_color_gradient2('Relative error (%)',      # titolo della legenda
                          low = 'red', mid = 'darkgreen', high = 'red',
                          midpoint = 0)+
    ylab("Relative error (%)")+
    xlab("wavelengths [nm]")
  
  save_path <- paste0(output_directory, filename,".png")
  
  ggsave(plot = plot, filename = save_path, height = 3.5, width = 7)
  
  return(plot)
}

#error_propagation ----
error_propagation <- function(std_hyp){
  length <- length(std_hyp)
  squares <- as.numeric(lapply(std_hyp,function(x){return(x^2)}))
  sum <- sum(squares)
  output <- sqrt(sum)/length
  return(output)
}

#error_propagation_v2 ----
error_propagation_v2 <- function(std_hyp, sample_std_hyp){
  length <- length(std_hyp)
  squares <- as.numeric(lapply(std_hyp,function(x){return(x^2)}))
  squares_with_sample_std <- as.numeric(lapply(squares,function(x){return(x + length*sample_std_hyp^2)})) #length riequilibra per ottenere sd quando non ho varianza sulle singole misure
  sum <- sum(squares_with_sample_std)
  output <- sqrt(sum)/length
  return(output)
}

#hypstar_spectral_resampling_ref_and_std ----
hypstar_spectral_resampling_ref_and_std <- function(sat_sensor,
                                                    dataset_res){
  #there is any need to consider propagation here?
  
  hypstar_std <- hypstar_spectral_resampling(sat_sensor = sat_sensor,
                                             dataset_res = dataset_res |> tidytable::select(-Ref) |> tidytable::rename(Ref = std)) |> 
    tidytable::mutate(sensor = sat_sensor) |>
    tidytable::rename(std = Ref)
  
  hypstar_Ref <- hypstar_spectral_resampling(sat_sensor = sat_sensor,
                                             dataset_res = dataset_res |> tidytable::select(-std)) |>
    tidytable::mutate(sensor = sat_sensor)
  
  hypstar_output <- hypstar_Ref |> tidytable::left_join(hypstar_std)
  
  return(hypstar_output)
}

#inverse_order_points ----
inverse_order_points <- function(x,y){
  exchange <- contoureR::orderPoints(cbind(x,y)[,1], cbind(x,y)[,2])
  seq <- 1:length(x)
  inverse <- data.frame(exchange,seq) |> tidytable::arrange(exchange) |> data.table::as.data.table()
  inverse_order <- inverse$seq
  return(inverse_order)
}

#save_function ----
save_function <- function(view_angle_output, filename, output_directory){
  write.table(x = view_angle_output,
              file = paste0(output_directory, filename),
              sep = ";", 
              quote = F, 
              row.names = F, 
              append = F)
}

#plot_polar ----
plot_polar <- function(time_or_sun_angle_output,
                       output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/",
                       title = paste0(unique(time_or_sun_angle_output$sensor),"_vs_hypstar",unique(time_or_sun_angle_output$operations)),
                       filename = paste0("polar_",gsub(";","_",title))){
  # library(tidytable)
  selected_data <- time_or_sun_angle_output |> tidytable::mutate(wvl_diff = abs(wvl - 900)) |>
    tidytable::filter(wvl_diff == min(wvl_diff)) |>
    tidytable::select(-wvl_diff)
  
  exact_wvl <- unique(selected_data$wvl)
  
  sat_coord <- selected_data |> tidytable::select(VAA_sat, VZA_sat) |> unique()
  
  library(ggplot2)
  plot <- ggplot(selected_data) +
    geom_point(aes(x = VAA_hyp, y = VZA_hyp, color = Ref_sat-Ref_hyp),
               size = 3) +
    coord_radial(start = -pi/2, end = +pi/2, inner.radius = 0.1,
                 r.axis.inside = T) +
    geom_point(data = sat_coord, aes(x = VAA_sat, y = VZA_sat),
               color = "black",
               size = 5,
               shape = "cross") +
    ylab("VZA [rad]") +
    xlab("VAA [rad]") +
    scale_color_gradient2('Ref_sat - Ref_hyp',      # titolo della legenda
                          low = 'blue', mid = 'green', high = 'red',
                          midpoint = 0)+
    ggtitle(paste0(title,"_at_",round(exact_wvl,0),"nm_","VZA_sat_",round(unique(selected_data$VZA_sat),0),"_VAA_sat_",round(unique(selected_data$VAA_sat),0)))
  
  save_path <- paste0(output_directory, filename,".png")
  
  ggsave(plot = plot, filename = save_path, height = 3.5, width = 10.5)
  
  return(plot)
}
