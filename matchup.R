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
  folder_list <- data.table::data.table(folder_name = list.dirs(home_path, recursive = FALSE, full.names = FALSE)) |>
    tidytable::mutate(date = as.Date(substr(folder_name, 4, 11), format = "%Y%m%d"),
           time = substr(folder_name, 12, 20)) |>
    tidytable::filter(date >= date_start,
           date <= date_end)
  
  hyp_function <- function(row_n, folder_list, home_path, bands_HYP){
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
    hyp_u_err_corr_systematic <- ncdf4::ncvar_get(ncin, "err_corr_systematic_reflectance")
    hyp_wvl  <- ncdf4::ncvar_get(ncin, "wavelength")
    hyp_vaa <-ncdf4::ncvar_get(ncin, "viewing_azimuth_angle")
    hyp_paa   <-ncdf4::ncvar_get(ncin, "pointing_azimuth_angle")
    hyp_vza  <-ncdf4::ncvar_get(ncin, "viewing_zenith_angle")
    hyp_sza <- ncdf4::ncvar_get(ncin, "solar_zenith_angle")
    hyp_saa <- ncdf4::ncvar_get(ncin, "solar_azimuth_angle")
    hyp_fwhm <- ncdf4::ncvar_get(ncin, "bandwidth")
    hyp_series_id <- ncdf4::ncvar_get(ncin, "series_id")
    
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
    
    complete_u_rel_random <- data.table::data.table(date = date,
                                                    time = time,
                                                    VZA = hyp_vza,
                                                    VAA = hyp_vaa,
                                                    PAA = hyp_paa,
                                                    SZA = hyp_sza,
                                                    SAA = hyp_saa) |>
      tidytable::bind_cols(hyp_u_rel_random |>
                             data.table::data.table() |>
                             `names<-`(paste0("B", 1:1554)) ) |>
      tidytable::mutate(across(c(VZA, VAA, PAA, SZA, SAA), ~round(.x)))
    
    complete_u_rel_random_long <- pivoting_longer_function(df = complete_u_rel_random,
                                                           bands = bands_HYP) |>
      tidytable::rename(u_rel_random = Ref)|>
      tidytable::mutate(u_rel_random = u_rel_random*0.01)
    
    complete_u_rel_systematic <- data.table::data.table(date = date,
                                                        time = time,
                                                        VZA = hyp_vza,
                                                        VAA = hyp_vaa,
                                                        PAA = hyp_paa,
                                                        SZA = hyp_sza,
                                                        SAA = hyp_saa) |>
      tidytable::bind_cols(hyp_u_rel_systematic |>
                             data.table::data.table() |>
                             `names<-`(paste0("B", 1:1554)) ) |>
      tidytable::mutate(across(c(VZA, VAA, PAA, SZA, SAA), ~round(.x)))
    
    complete_u_rel_systematic_long <- pivoting_longer_function(df = complete_u_rel_systematic,
                                                               bands = bands_HYP) |>
      tidytable::rename(u_rel_systematic = Ref) |>
      tidytable::mutate(u_rel_systematic = u_rel_systematic*0.01) #rescale factor
    
    # complete_err_corr_systematic <- data.table::data.table(date = date,
    #                                                     time = time,
    #                                                     VZA = hyp_vza,
    #                                                     VAA = hyp_vaa,
    #                                                     PAA = hyp_paa,
    #                                                     SZA = hyp_sza,
    #                                                     SAA = hyp_saa) |>
    #   tidytable::bind_cols(hyp_u_err_corr_systematic |>
    #                          data.table::data.table() |>
    #                          `names<-`(paste0("B", 1:1554)) ) |>
    #   tidytable::mutate(across(c(VZA, VAA, PAA, SZA, SAA), ~round(.x))) 
    # 
    # complete_err_corr_systematic_long <- pivoting_longer_function(df = complete_err_corr_systematic,
    #                                                            bands = bands_HYP) |>
    #   tidytable::rename(err_corr_systematic = Ref) |>
    #   tidytable::mutate(err_corr_systematic = err_corr_systematic*0.01) #rescale factor
    
    complete_DB <- complete_ref_long |> tidytable::left_join(complete_std_long) |>
      tidytable::left_join(complete_u_rel_random_long) |> tidytable::left_join(complete_u_rel_systematic_long) |>
      #tidytable::left_join(complete_err_corr_systematic_long) |>
      tidytable::mutate(u_rand = u_rel_random*Ref,
                        u_sys = u_rel_systematic*Ref) |>
      tidytable::select(-c(u_rel_random,u_rel_systematic))
    
    return(complete_DB)
    
  }
  
  # print("start_reading")
  
  res <- mirai::mirai_map(1:nrow(folder_list), function(x){return(hyp_function(x,folder_list,home_path,bands_HYP))}, 
                          hyp_function = hyp_function, 
                          folder_list = folder_list, 
                          home_path = home_path, 
                          bands_HYP = bands_HYP,
                          pivoting_longer_function = pivoting_longer_function)[] |> tidytable::bind_rows() 
  
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
  output <- filtered  |> tidytable::select(-PAA, -u_rand, -u_sys) |> #al momento li tolgo 
    tidytable::mutate(time = substring(time,2,8)) |>
    tidytable::mutate(time = sub("\\s+$", "", gsub('(.{2})', '\\1:', time))) |>
    tidytable::mutate(time = hms::as_hms(substring(time,1,8))) |>
    tidytable::mutate(hyp_sensor = "JSIT")
  
  if(is.null(filename)){
    save_function(output, paste0("hypstar_data_from_", date_start, "_to_", date_end, ".csv"), output_directory)
  }else{
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

#pixel_aggregation_function ----
pixel_aggregation_function <- function(function_type = "exact", #choose between "exact" and "terra"
                                       aggregation_type = "mean", #choose between "mean" and "weight"
                                       raster, 
                                       vector,
                                       max_number_of_pixels){
  if(function_type == "exact"){
    #pixel extraction done using exact_extract
    pixel_extraction <- exactextractr::exact_extract(raster, vector, coverage_area = F,
                                                     include_cols = names(vector)[names(vector) != "geometry"]) |>
      tidytable::bind_rows() |>
      tidytable::select(-FID,-X,-y)
  }else if(function_type == "terra"){
    #pixel extraction done using terra
    pixel_extraction <- terra::extract(raster, vector, fun = "table", touches = T, weights = T) |>
      tidytable::select(-ID) |>
      `names<-`(c("value","coverage_fraction"))
  }else{
    stop("Error, choose a function_type among 'exact' and 'terra'")
  }
  
  if(aggregation_type == "mean"){
    #here we perform the simple mean using all the max_number_of_pixels closest pixels  
    pixel_aggregation <- pixel_extraction |>
      tidytable::arrange(desc(coverage_fraction)) |>
      tidytable::select(-coverage_fraction) |>
      tidytable::filter(tidytable::row_number() <= max_number_of_pixels) |>
      tidytable::summarise(sd = sd(value, na.rm = T), mean = mean(value, na.rm=T))
  }else if(aggregation_type == "weight"){
    #here we perform the weighted mean using all the max_number_of_pixels closest pixels
    pixel_aggregation <- pixel_extraction |>
      tidytable::arrange(desc(coverage_fraction)) |>
      tidytable::filter(tidytable::row_number() <= max_number_of_pixels) |>
      tidytable::summarise(mean = weighted.mean(x = value, w = coverage_fraction, na.rm=T),
                           sd = sd(value, na.rm = T))
  }else{
    stop("Error, choose an aggregation_type among 'mean' and 'weight'.")
  }
  
  return(pixel_aggregation)
}

#get_sat_date ----
get_sat_date <- function(satellite,product_path){
  if(satellite == "S2"){
    return(as.Date(substring(basename(product_path),8,15),"%Y%m%d"))
  }
    
  if(satellite == "PRS"){
    read_angle_and_time <- read.csv(list.files(path = product_path,
                                               pattern = "all_angles_file.csv",
                                               full.names = T,
                                               recursive = T))
    return(as.Date(substring(read_angle_and_time$date, 1,10),"%Y-%m-%d"))
  }
  if(satellite == "En"){
    return(as.Date(substring(basename(product_path),30,37),"%Y%m%d"))
  }
}

#get_sat_time ----
get_sat_time <- function(satellite,product_path){
  if(satellite == "S2"){
    return(hms::as_hms(substring(sub("\\s+$", "", gsub('(.{2})', '\\1:',substring(basename(product_path),17,22))),1,8)))
  }
  
  if(satellite == "PRS"){
    read_angle_and_time <- read.csv(list.files(path = product_path,
                                               pattern = "all_angles_file.csv",
                                               full.names = T,
                                               recursive = T))
    return(hms::as_hms(substring(read_angle_and_time$date, 12,19)))
  }
  
  if(satellite == "En"){
    return(hms::as_hms(substring(sub("\\s+$", "", gsub('(.{2})', '\\1:', substring(basename(product_path),39,44))),1,8)))
    
  }
    
}

#get_sat_sensor_ID ----
get_sat_sensor_ID <- function(satellite,product_path){
  if(satellite == "S2"){
    return(substring(basename(product_path),24,26))
  }
  
  if(satellite == "PRS"){
    return("PRISMA")
  }
  
  if(satellite == "En"){
    return("EnMAP")
  }
  
  stop("Error, not valid satellite")
}

#get_sat_bands_list_and_names ----
get_sat_bands_list_and_names <- function(satellite,product_path){
  if(satellite == "S2"){
    list_of_tifs <- list.files(product_path, full.names = T, recursive = T,
                               pattern = "\\.tif$")
    
    bands <- list_of_tifs[!stringr::str_detect(list_of_tifs, "mask")]
    
    list_of_rasters <- lapply(bands, terra::rast)
    
    extract_names <- stringr::str_match(bands, "m_\\s*(.*?)\\s*.tif")[,2]
    
    return(list(list_of_rasters,extract_names))
  }
  
  if(satellite == "PRS"){
    if(grepl(pattern = "L1", product_path)){
      search_for = "CARST"
    }else if(grepl(pattern = "L2", product_path)){
      search_for = "CRST"
    }
    image_path <- list.files(product_path, full.names = T, recursive = T, pattern = paste0("\\",search_for,".tif$"))
    
    print(image_path)
    
    image <- terra::rast(image_path)
    
    list_of_rasters <- lapply(1:terra::nlyr(image), function(x){return(image[[x]])})
    
    extract_names <- unlist(lapply(strsplit(names(image)," "), function(x){return(x[[1]])}))
    
    return(list(list_of_rasters,extract_names))
  }
  
  if(satellite == "En"){
    print(product_path)
    
    image_path <- list.files(paste0(product_path, "/"), full.names = T, recursive = T,
                             pattern = "\\.BSQ$")
    
    image <- terra::rast(image_path)
    
    list_of_rasters <- lapply(1:terra::nlyr(image), function(x){return(image[[x]])})
    
    extract_names <- unlist(lapply(strsplit(names(image)," "), function(x){return(paste0("B",x[[2]]))}))
    return(list(list_of_rasters,extract_names))
  }
}

#get_sat_angles ----
get_sat_angles <- function(satellite, product_path, bands_statistics, bandnames){
  if(satellite == "S2"){
    angles <- bands_statistics |>
      tidytable::select(mean) |>
      tidytable::rename(value = mean) |>
      tidytable::mutate(bands = bandnames, 
                        value = as.numeric(value)) |>
      tidytable::filter(!grep(bands, pattern = "^B.*")) |>
      tidytable::pivot_wider(names_from = bands, values_from = value)
  }
  if(satellite == "PRS"){
    read_angle_and_time <- tidytable::fread(list.files(path = product_path,
                                                       pattern = "all_angles_file.csv",
                                                       full.names = T,
                                                       recursive = T))
    
    angles <- read_angle_and_time |>
      tidytable::rename(SAA = sunaz,
                        SZA = sunzen,
                        VZA = sensor_zenith,
                        VAA = sensor_azimuth) |>
      tidytable::select(VAA,VZA,SZA,SAA)
  }
  if(satellite == "En"){
    metadata_path <- list.files(paste0(product_path, "/"), full.names = T, recursive = T,
                                pattern = "\\METADATA.XML$")
    
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
    
    angles <- data.frame(SAA = sunAzimithAngle,
                         SZA = sunZenithAngle,
                         VZA = viewingZenithAngle,
                         VAA = viewingAzimithAngle) |>
      tidytable::as_tidytable()
  }
  return(angles)
}

#get_statistics_dataset ----
get_statistics_dataset <- function(satellite,
                                   product_path, 
                                   hypstar_point_buffered, 
                                   aggregation_type, 
                                   number_closest_pixels, 
                                   function_type){
  
  print(product_path)
  
  bands_data <- get_sat_bands_list_and_names(satellite = satellite,
                                             product_path = product_path)
  
  date_ID <- get_sat_date(satellite = satellite,
                       product_path = product_path)
  
  time_ID <- get_sat_time(satellite = satellite,
                       product_path = product_path)
  
  sensor_ID <- get_sat_sensor_ID(satellite = satellite,
                                 product_path = product_path)
  
  bands_statistics <- lapply(bands_data[[1]], pixel_aggregation_function, vector = hypstar_point_buffered,
                             max_number_of_pixels = number_closest_pixels, aggregation_type = aggregation_type,
                             function_type = function_type) |>
    tidytable::bind_rows()
  
  angles <- get_sat_angles(satellite = satellite, product_path = product_path, bands_statistics = bands_statistics,
                           bandnames = bands_data[[2]])

  statistics_pivoted <- bands_statistics |> 
    tidytable::mutate(bands = bands_data[[2]]) |>
    tidytable::rename(Ref = mean) |>
    tidytable::filter(grep(bands, pattern = "^B.*")) |>
    tidytable::mutate(date = date_ID) |>
    tidytable::mutate(time = time_ID) |>
    tidytable::mutate(sensor = sensor_ID) |>
    tidytable::bind_cols(angles) |>
    data.table::as.data.table()
  
  return(statistics_pivoted)
}

#get_s2_data ----
get_s2_data <- function(date_start, 
                        date_end,
                        hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                        repo_root_path = "//10.0.1.243/Repository/nr_products/sentinel2/l2/T32TQQ/", #masked S2 data
                        filename = NULL,
                        number_closest_pixels = 36,
                        output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/",
                        scale = F,
                        aggregation_type = "mean", # this can be 1,2,3 depending on: 
                        function_type = "exact"){
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
  
  print(product_paths)
  
  hypstar_point <- sf::st_read(hypstar_point_path)
  
  hypstar_point_buffered <- sf::st_buffer(hypstar_point, dist = 30) 
  
  all_datasets <- mirai::mirai_map(product_paths, function(x){return(get_statistics_dataset(satellite,x,hypstar_point_buffered,aggregation_type,number_closest_pixels,function_type))},
                                   get_statistics_dataset = get_statistics_dataset,
                                   pixel_aggregation_function = pixel_aggregation_function,
                                   hypstar_point_buffered = hypstar_point_buffered,
                                   aggregation_type = aggregation_type, 
                                   number_closest_pixels = number_closest_pixels, 
                                   function_type = function_type,
                                   satellite = "S2",
                                   get_sat_date = get_sat_date,
                                   get_sat_time = get_sat_time,
                                   get_sat_sensor_ID = get_sat_sensor_ID,
                                   get_sat_bands_list_and_names =get_sat_bands_list_and_names,
                                   get_sat_angles = get_sat_angles)[]
  
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
  }else{
    save_function(with_wvl, paste0("S2_data_from_",date_start,"_to_",date_end, ".csv"), output_directory)
  }
  
  return(with_wvl)
}


#get_PRISMA_data ----
get_PRISMA_data <- function(date_start, 
                            date_end,
                            hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                            repo_root_path = "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/processing_L1", # o //10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/put_PRISMA_he5_and_S2_tif_here_L2
                            filename = NULL,
                            output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/",
                            number_closest_pixels = 4,
                            scale = F,
                            aggregation_type = "mean", 
                            function_type = "exact"){
  year_start <- lubridate::year(as.Date(date_start))
  year_end <- lubridate::year(as.Date(date_end))
  years <- as.character(year_start:year_end)
  repo_year_path <- paste0(repo_root_path, years, "/", PRISMA_level_product, "/")
  
  #now you enter in the folder where there are products
  #TEMPORARY
  list_of_products <- list.dirs(path = repo_root_path,full.names = T,recursive = F)
  list_of_dates <- lapply(list_of_products,function(x){return(as.Date(substring(strsplit(basename(x), split = "_")[[1]][[5]],"1","8"), format = "%Y%m%d"))})
  product_paths <- unlist(list_of_products[list_of_dates > as.Date(date_start) & list_of_dates < as.Date(date_end)])
  
  hypstar_point <- sf::st_read(hypstar_point_path)
  
  hypstar_point_buffered <- sf::st_buffer(hypstar_point, dist = 30) 
  
  all_datasets <- mirai::mirai_map(product_paths, function(x){return(get_statistics_dataset(satellite,x,hypstar_point_buffered,aggregation_type,number_closest_pixels,function_type))},
                                   get_statistics_dataset = get_statistics_dataset,
                                   pixel_aggregation_function = pixel_aggregation_function,
                                   hypstar_point_buffered = hypstar_point_buffered,
                                   aggregation_type = aggregation_type, 
                                   number_closest_pixels = number_closest_pixels, 
                                   function_type = function_type,
                                   satellite = "PRS",
                                   get_sat_date = get_sat_date,
                                   get_sat_time = get_sat_time,
                                   get_sat_sensor_ID = get_sat_sensor_ID,
                                   get_sat_bands_list_and_names =get_sat_bands_list_and_names,
                                   get_sat_angles = get_sat_angles)[]
  
  binded <- dplyr::bind_rows(all_datasets)|>
    stats::na.omit()
  
  if(scale){
    binded <- binded |> 
      tidytable::mutate(Ref = Ref/10000) |>
      tidytable::mutate(sd = sd/10000)
  }
  
  with_wvl <- add_wvl(binded) |>
    stats::na.omit() 
  
  if(is.null(filename)){
    save_function(with_wvl, paste0(paste0(basename(product_paths), collapse = "_and_"), ".csv"), output_directory)
  }else{
    save_function(with_wvl, filename, output_directory)
  }
  
  return(with_wvl)
  
}

#get_EnMAP_data ----
get_EnMAP_data <- function(date_start, 
                           date_end,
                           EnMAP_level_product = "L2A", #enmap level product
                           hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp", #hypstar shapefile path
                           repo_root_path = "//10.0.1.243/nr_data/3_rs_data/EnMAP/JDS/", #the main folder of EnMAP data
                           number_closest_pixels = 4,
                           filename = NULL, #if not NULL the enmap data will be saved with this filename
                           output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/", #output directory if filename is not null
                           scale = F,
                           aggregation_type = "mean", 
                           function_type = "exact"){
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
  
  hypstar_point_buffered <- sf::st_buffer(hypstar_point, dist = 30) 
  
  all_datasets <- mirai::mirai_map(product_paths, function(x){return(get_statistics_dataset(satellite,x,hypstar_point_buffered,aggregation_type,number_closest_pixels,function_type))},
                                   get_statistics_dataset = get_statistics_dataset,
                                   pixel_aggregation_function = pixel_aggregation_function,
                                   hypstar_point_buffered = hypstar_point_buffered,
                                   aggregation_type = aggregation_type, 
                                   number_closest_pixels = number_closest_pixels, 
                                   function_type = function_type,
                                   satellite = "En",
                                   get_sat_date = get_sat_date,
                                   get_sat_time = get_sat_time,
                                   get_sat_sensor_ID = get_sat_sensor_ID,
                                   get_sat_bands_list_and_names =get_sat_bands_list_and_names,
                                   get_sat_angles = get_sat_angles)[]
  
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
  }else{
    save_function(with_wvl, paste0(paste0(basename(product_paths), collapse = "_and_"), ".csv"), output_directory)
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
    plot_divide_by_date(plot_polar, dataset = time_or_sun_angle_output)
    # plot_polar(time_or_sun_angle_output = time_or_sun_angle_output,
    #            output_directory = gsub("csv","plots", output_directory))
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
    save_function(view_angle_output, paste0(sat_sensor,"_matchups_from_",min(sat_data$date),"_to_",max(sat_data$date),
                                            "_TIME-or-ANGLE-",
                                            strsplit(strsplit(unique(view_angle_output$operations), ":")[[1]][[2]],"_")[[1]][[1]],
                                            "_WINDOW-",
                                            invisible(suppressWarnings(as.numeric(tryCatch(readr::parse_number(strsplit(unique(view_angle_output$operations), ":")[[1]][[2]]),
                                                                error = function(e){return(NULL)})))),
                                            strsplit(strsplit(unique(view_angle_output$operations), ":")[[1]][[2]],"_")[[1]][[length(strsplit(strsplit(unique(view_angle_output$operations), ":")[[1]][[2]],"_")[[1]])]],
                                            "_VIEW-",
                                            strsplit(strsplit(unique(view_angle_output$operations), ":")[[1]][[3]],"_")[[1]][[1]],
                                            "_WINDOW-",
                                            invisible(suppressWarnings(as.numeric(tryCatch(readr::parse_number(strsplit(unique(view_angle_output$operations), ":")[[1]][[3]]),
                                                                error = function(e){return(NULL)})))),
                                            "deg",
                                            ".csv"), output_directory) 
  }
  
  
  return(view_angle_output)
  
}

#matchup_function_v2 ----
matchup_function_v2 <- function(sat_data, 
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
  
  
  #DATE PROBLEM
  date_output <- date_function(sat_col_renamed = sat_col_renamed,
                               hypstar_col_renamed = hypstar_col_renamed,
                               only_same_date = only_same_date)
  
  #GEOMETRY PROBLEM
  geometry_output <- elemental_matchup(date_output = date_output)
  
  #save
  if(!is.null(filename)){
    save_function(geometry_output, filename, output_directory)
  }else{
    save_function(geometry_output, paste0(sat_sensor,"_matchups_from_",min(sat_data$date),"_to_",max(sat_data$date),
                                            "_elemental",
                                            ".csv"), output_directory) 
  }
  
  return(geometry_output)
}

#elemental_matchup ----
elemental_matchup <- function(date_output){
  compute_best_config <- date_output |> tidytable::mutate(rel_azimuth_sat = abs(SAA_sat-VAA_sat),
                                                                rel_azimuth_hyp = abs(SAA_hyp-VAA_hyp)) |>
    #using now the cosine law https://en.wikipedia.org/wiki/Law_of_cosines
    # tidytable::mutate(rel_angle_sat = acos( (cos(a*pi/180)/cos(b*pi/180) + cos(b*pi/180)/cos(a*pi/180) - tan(a*pi/180)*sin(a*pi/180)*cos(b*pi/180) - tan(b*pi/180)*sin(b*pi/180)*cos(a*pi/180) + 2*sin(a*pi/180)*sin(b*pi/180)*cos(c*pi/180) )/2) *180/pi ) |>
    tidytable::mutate(rel_angle_sat = acos( (cos(VZA_sat*pi/180)/cos(SZA_sat*pi/180) + cos(SZA_sat*pi/180)/cos(VZA_sat*pi/180) - tan(VZA_sat*pi/180)*sin(VZA_sat*pi/180)*cos(SZA_sat*pi/180) - tan(SZA_sat*pi/180)*sin(SZA_sat*pi/180)*cos(VZA_sat*pi/180) + 2*sin(VZA_sat*pi/180)*sin(SZA_sat*pi/180)*cos(rel_azimuth_sat*pi/180) )/2) *180/pi ) |>
    tidytable::mutate(rel_angle_hyp = acos( (cos(VZA_hyp*pi/180)/cos(SZA_hyp*pi/180) + cos(SZA_hyp*pi/180)/cos(VZA_hyp*pi/180) - tan(VZA_hyp*pi/180)*sin(VZA_hyp*pi/180)*cos(SZA_hyp*pi/180) - tan(SZA_hyp*pi/180)*sin(SZA_hyp*pi/180)*cos(VZA_hyp*pi/180) + 2*sin(VZA_hyp*pi/180)*sin(SZA_hyp*pi/180)*cos(rel_azimuth_hyp*pi/180) )/2) *180/pi ) |>
    #using now the dihedral angle https://en.wikipedia.org/wiki/Dihedral_angle
    # tidytable::mutate(plane_inclination = acos(abs(tan(a*pi/180))/sqrt(1+tan(a*pi/180)^2+((tan(a*pi/180)-tan(b*pi/180)*cos(c*pi/180))/(tan(b*pi/180)*sin(c*pi/180)))^2))*180/pi ) |>
    tidytable::mutate(plane_inclination_sat = acos(abs(tan(VZA_sat*pi/180))/sqrt(1+tan(VZA_sat*pi/180)^2+((tan(VZA_sat*pi/180)-tan(SZA_sat*pi/180)*cos(rel_angle_sat*pi/180))/(tan(SZA_sat*pi/180)*sin(rel_angle_sat*pi/180)))^2))*180/pi ) |>
    tidytable::mutate(plane_inclination_hyp = acos(abs(tan(VZA_hyp*pi/180))/sqrt(1+tan(VZA_hyp*pi/180)^2+((tan(VZA_hyp*pi/180)-tan(SZA_hyp*pi/180)*cos(rel_angle_hyp*pi/180))/(tan(SZA_hyp*pi/180)*sin(rel_angle_hyp*pi/180)))^2))*180/pi ) |>
    tidytable::group_by(date, time_hyp, VZA_hyp, VAA_hyp) |>
    tidytable::mutate(plane_inclination_diff = abs(plane_inclination_sat-plane_inclination_hyp)) |>
    tidytable::mutate(rel_angle_diff = abs(rel_angle_sat-rel_angle_hyp)) |>
    tidytable::ungroup() |>
    tidytable::filter(rel_angle_diff == min(rel_angle_diff)) |>
    tidytable::select(-rel_angle_diff,-rel_azimuth_sat,-rel_azimuth_hyp,-rel_angle_sat,-rel_angle_hyp,
                      -plane_inclination_sat,-plane_inclination_hyp,
                      -plane_inclination_diff) |>
    tidytable::mutate(operations = ":elemental:elemental:elemental:elemental") 
  
  # tidytable::fwrite(compute_best_config,
  #                   "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/prova_configurazione_elementale/configurazione_elementale.xlsx")
  
  #la divisione in rel_angle_diff e plane_inclination è motivata dalla seguente considerazione
  #
  #rel_angle_diff riguarda (principalmente) la microstruttura: 
  #in base alla geometria di diffusione della singola particella, guardare con lo stesso angolo 
  #relativo mostra lo stesso spettro di riflessione. 
  #Non è possibile escludere del tutto la macrostruttura a causa delle riflessioni multiple nel 
  #target
  #
  #plane_inclination riguarda solamente la macrostruttura:
  #cambiare piano di inclinazione mostra (fa una radiografia) della diversa macrostruttura
  #del target, ovvero di come le particelle si aggregano nel target. La singola particella infatti non 
  #può distinguere il piano di inclinazione dei raggi pertanto la microstruttura cambiando 
  #il piano di inclinazione dei raggi non può essere analizzata
  
  return(compute_best_config)
}


#view_mean ----
view_mean <- function(med_view_angle_output, range_view_angle, weighted_mean){
  cols <- c(operations = "")
  cols_num <- c(number_data_in_window = "")
  
  write_ops <- paste0("average_in_view_angle_window_of",range_view_angle,"_degrees")
  
  if(weighted_mean){
    write_ops <- paste0("weighted",write_ops)
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
    tidytable::mutate(number_data_in_window = paste(number_data_in_window,length(std_hyp),sep=":")) |>
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
      tidytable::mutate(operations = paste(operations,write_ops,sep=":") )
    
  }else{
    #as requested from GC when we do the weighted mean the viewing angle geometry recorded is the one
    #of the satellite
    
    names_cols <- names(ref_averaged)
    
    view_angle_output <- ref_averaged |>
      tibble::add_column(!!!cols[!names(cols) %in% names_cols]) |>
      tidytable::mutate(operations = paste(operations,write_ops,sep=":") ) |> 
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
    tidytable::mutate(operations = paste(operations,write_ops,sep=":") ) 
  
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
    tidytable::mutate(operations = paste(operations,write_ops,sep=":") ) 
  
  return(closest_output_check)
}

#time_or_sun_mean ----
time_or_sun_mean <- function(med_output, range_time_or_sun_angle, time_or_sun_angle, weighted_mean){
  
  cols <- c(operations = "")
  cols_num <- c(number_data_in_window = "")
  
  #filter over time/angle (keep only times/angles in the range)
  if(time_or_sun_angle){
    write_ops <- paste0("average_in_time_window_of_",range_time_or_sun_angle,"_hours")
    filtered <- med_output |>
      tidytable::filter(time_diff < range_time_or_sun_angle/2*3600) |>
      tidytable::mutate(weight = 1/as.numeric(abs(time_diff))) |>
      tidytable::select(-time_diff) 
  }else{
    write_ops <- paste0("average_in_sun_angle_window_of",range_time_or_sun_angle,"_degrees")
    filtered <- med_output |>
      tidytable::filter(angle_diff < range_time_or_sun_angle/2) |>
      tidytable::mutate(weight = 1/as.numeric(abs(angle_diff))) |>
      tidytable::select(-angle_diff)
  }
  
  if(weighted_mean){
    write_ops <- paste0("weighted",write_ops)
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
    tidytable::mutate(number_data_in_window = paste(number_data_in_window,length(std_hyp),sep=":")) |>
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
      tidytable::mutate(operations = paste(operations,write_ops,sep=":") ) |> 
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
      tidytable::mutate(operations = paste(operations,write_ops,sep=":") ) |> 
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
                        output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/",
                        filename = NULL,
                        limits = NULL){
  # library(ggplot2)
  
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
  
  plot_function <- function(date_sel, limits, matchup_result, output_directory, filename, sat_color, alpha,color_station){
    data_plot <- matchup_result |> tidytable::filter(date == date_sel)
    title <- make_title_for_matchup_plots(data_plot)
    plot <- ggplot(data = data_plot) + 
      geom_line(aes(x = wvl, y = Ref_sat), linewidth = 1, color = sat_color) +
      geom_line(aes(x = wvl, y = Ref_hyp), linewidth = 1, color = color_station) +
      geom_ribbon(aes(x = wvl, y = Ref_sat, ymin = Ref_sat - sd_sat, ymax = Ref_sat + sd_sat), alpha = alpha, fill = sat_color) +
      geom_ribbon(aes(x = wvl, y = Ref_hyp, ymin = Ref_hyp - std_hyp, ymax = Ref_hyp + std_hyp), alpha = alpha, fill = color_station) +
      ggtitle(title) +
      ylab("R [-]") +
      xlab("Wavelength [nm]") 
    
    if(!is.null(limits)){
      plot <- plot + ylim(limits)
    }
    
    if(is.null(filename)){
      filename <- paste0("ribbon_",title)
    }
    
    save_path <- paste0(output_directory, filename,".png")
    
    ggsave(plot = plot, filename = save_path, height = 3.5, width = 14)
    
    return(plot)
  }
  
  mirai::everywhere(
    library(ggplot2)
  )
  
  val <- mirai::mirai_map(unique(matchup_result$date),
                   function(x){return(plot_function(x,
                                                    limits,
                                                    matchup_result,
                                                    output_directory,
                                                    filename,
                                                    sat_color,
                                                    alpha,
                                                    color_station))},
                   plot_function = plot_function,
                   limits = limits,
                   matchup_result = matchup_result,
                   output_directory = output_directory,
                   filename = filename,
                   sat_color = sat_color,
                   alpha = alpha,
                   color_station = color_station,
                   make_title_for_matchup_plots = make_title_for_matchup_plots)[]
  
  return(val)
}

#plot_rainbow ----
plot_rainbow <- function(matchup_result,
                         point_size = 3,
                         output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/",
                         filename = NULL){
  # library(ggplot2)
  palette <- c(rev(rainbow(8)[-length(rainbow(8))]),"#bbbbbb","#777777","#444444","#000000")
  values_palette <- scales::rescale(x = c(350,400,450,500,550,600,650,700,1000,1300,1800),
                            to = c(0,1))
  sat_name <- unique(substr(matchup_result$sensor,1,2))
  
  plot_function <- function(data_sel, matchup_result, point_size, sat_name, values_palette, palette, output_directory, filename){
    data_plot <- matchup_result |> tidytable::filter(date == data_sel)
    title <- make_title_for_matchup_plots(data_plot)
    plot <- ggplot() +
      geom_point(data = data_plot, 
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
    
    if(is.null(filename)){
      filename <- paste0("rainbow_",title)
    }
    save_path <- paste0(output_directory, filename,".png")
    
    ggsave(plot = plot, filename = save_path, height = 5, width = 13)
    
    return(plot)
  }
  
  mirai::everywhere(
    library(ggplot2)
  )
  val <- mirai::mirai_map(unique(matchup_result$date),
                          function(x){return(plot_function(x, matchup_result, point_size, sat_name, values_palette, palette, output_directory, filename))},
                          plot_function = plot_function,
                          matchup_result = matchup_result, 
                          point_size = point_size, 
                          sat_name = sat_name, 
                          values_palette = values_palette, 
                          palette = palette,
                          output_directory = output_directory, 
                          filename = filename,
                          make_title_for_matchup_plots = make_title_for_matchup_plots)[]
  
  
  return(val)
}

#plot_relative_errors ----
plot_relative_errors <- function(matchup_result,
                                 point_size = 3,
                                 
                                 output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/",
                                 filename = NULL){
  relative_error_df <- matchup_result |> tidytable::mutate(relative_error = (Ref_sat-Ref_hyp)/ Ref_hyp *100) 
  # library(ggplot2)
  
  plot_function <- function(date_sel, point_size, output_directory, filename , relative_error_df){
    
    data_plot <- relative_error_df |> tidytable::filter(date == date_sel)
    title <- make_title_for_matchup_plots(data_plot)
    
    plot <- ggplot(data = data_plot) +
      geom_abline(slope = 0, 
                  intercept = 0,
                  linetype = "dashed") +
      ggtitle(title)+
      geom_point(mapping = aes(x = wvl, y = relative_error, color = relative_error), size = point_size)+
      scale_color_gradient2('Relative error (%)',      # titolo della legenda
                            low = 'red', mid = 'darkgreen', high = 'red',
                            midpoint = 0)+
      ylab("Relative error (%)")+
      xlab("wavelengths [nm]")
    
    if(is.null(filename)){
      filename <- paste0("relative_error_",title)
    }
    
    save_path <- paste0(output_directory, filename,".png")
    
    ggsave(plot = plot, filename = save_path, height = 5, width = 13)
    
    return(plot)
  }
  mirai::everywhere(
    library(ggplot2)
  )
  val <- mirai::mirai_map(unique(relative_error_df$date),
                          function(x){return(plot_function(x, point_size, output_directory, filename, relative_error_df))},
                          plot_function = plot_function,
                          make_title_for_matchup_plots = make_title_for_matchup_plots,
                          point_size = point_size, 
                          output_directory = output_directory, 
                          filename = filename, 
                          relative_error_df = relative_error_df)[]
  
  
  return(val)
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
  squares_with_sample_std <- as.numeric(lapply(std_hyp,function(x){return(x^2 + length*sample_std_hyp^2)})) #length riequilibra per ottenere sd quando non ho varianza sulle singole misure
  return(sqrt(sum(squares_with_sample_std))/length)
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
                       filename = NULL){
  selected_data <- time_or_sun_angle_output |> tidytable::mutate(wvl_diff = abs(wvl - 900)) |>
      tidytable::filter(wvl_diff == min(wvl_diff)) |>
      tidytable::select(-wvl_diff) 
    
    exact_wvl <- unique(selected_data$wvl)
    
    sat_coord <- selected_data |> tidytable::select(VAA_sat, VZA_sat) |> unique()
    
    library(ggplot2)
    plot <- ggplot(selected_data) +
      geom_point(aes(x = VAA_hyp, y = VZA_hyp, color = Ref_sat-Ref_hyp),
                 size = 3) +
      coord_radial(theta = "x", start = -pi*15/180, thetalim = c(0,330), inner.radius = 0.1,
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
      scale_x_continuous(
        limits = c(0, 360),
        breaks = seq(0, 330, by = 30),
        minor_breaks = NULL
      )+
      ggtitle(label = "polar_plot", subtitle = paste0(make_title_for_polar_plots(selected_data),"_at_",round(exact_wvl,0),"nm_","VZA_sat_",round(unique(selected_data$VZA_sat),0),"_VAA_sat_",round(unique(selected_data$VAA_sat),0)))
    
    if(is.null(filename)){
      filename <- paste0("polar_plot_",paste0(make_title_for_polar_plots(selected_data),"_at_",round(exact_wvl,0),"nm_","VZA_sat_",round(unique(selected_data$VZA_sat),0),"_VAA_sat_",round(unique(selected_data$VAA_sat),0)))
    }
    
    save_path <- paste0(output_directory, filename,".png")
    
    ggsave(plot = plot, filename = save_path, height = 9, width = 19)

  return(plot)
}

#plot_all_view_spectra ----
plot_all_view_spectra <- function(time_or_sun_angle_output,
                                  output_directory = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/plots/",
                                  title = paste0(unique(time_or_sun_angle_output$sensor),"_vs_hypstar",unique(time_or_sun_angle_output$operations)),
                                  filename = paste0("all_view_spectra_",gsub(":","_",title))){
  
}

#make_title_for_polar_plots ----
make_title_for_polar_plots <- function(time_or_sun_angle_output){
  title = paste0("SAT-",unique(substr(time_or_sun_angle_output$sensor,1,2)),
                 "_HYP-",
                 unique(substr(time_or_sun_angle_output$hyp_sensor,1,4)),
                 "_DT-", 
                 gsub("-","",unique(time_or_sun_angle_output$date)),
                 "T",
                 gsub(":","",unique(time_or_sun_angle_output$time_sat)),
                 "_TIME-or-ANGLE-",
                 strsplit(strsplit(unique(time_or_sun_angle_output$operations), ":")[[1]][[2]],"_")[[1]][[1]],
                 "_WINDOW-",
                 invisible(suppressWarnings(as.numeric(tryCatch(readr::parse_number(strsplit(unique(time_or_sun_angle_output$operations), ":")[[1]][[2]]),
                                                                error = function(e){return(NULL)})))),
                 strsplit(strsplit(unique(time_or_sun_angle_output$operations), ":")[[1]][[2]],"_")[[1]][[length(strsplit(strsplit(unique(time_or_sun_angle_output$operations), ":")[[1]][[2]],"_")[[1]])]]
  )
  return(title)
}

#make_title_for_matchup_plots ----
make_title_for_matchup_plots <- function(matchup_result){
  title = paste0("SAT-",unique(substr(matchup_result$sensor,1,2)),
                 "_HYP-",
                 unique(substr(matchup_result$hyp_sensor,1,4)),
                 "_DT-", 
                 gsub("-","",unique(matchup_result$date)),
                 "T",
                 gsub(":","",unique(matchup_result$time_sat)),
                 "_TIME-or-ANGLE-",
                 strsplit(strsplit(unique(matchup_result$operations), ":")[[1]][[2]],"_")[[1]][[1]],
                 "_WINDOW-",
                 invisible(suppressWarnings(as.numeric(tryCatch(readr::parse_number(strsplit(unique(matchup_result$operations), ":")[[1]][[2]]),
                          error = function(e){return(NULL)})))),
                 strsplit(strsplit(unique(matchup_result$operations), ":")[[1]][[2]],"_")[[1]][[length(strsplit(strsplit(unique(matchup_result$operations), ":")[[1]][[2]],"_")[[1]])]],
                 "_VIEW-",
                 strsplit(strsplit(unique(matchup_result$operations), ":")[[1]][[3]],"_")[[1]][[1]],
                 "_WINDOW-",
                 invisible(suppressWarnings(as.numeric(tryCatch(readr::parse_number(strsplit(unique(matchup_result$operations), ":")[[1]][[3]]),
                                     error = function(e){return(NULL)})))),
                 "deg"
                 )
  return(title)
}

#get_list_of_date_datasets ----
get_list_of_date_datasets <- function(date_to_keep, dataset){
  return(dataset |> tidytable::filter(date == date_to_keep))
}

#plot_divide_by_date ----
plot_divide_by_date <- function(fns_to_plot, dataset){
  list_of_dates = unique(dataset$date)
  
  list_of_date_datasets <- lapply(list_of_dates, get_list_of_date_datasets, dataset = dataset)
  
  lapply(list_of_date_datasets, fns_to_plot)
}
