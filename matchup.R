#Create match-up
#get_hypstar_data ----
get_hypstar_data <- function(date_start,
                             date_end,
                             home_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/PROCESSED",
                             output_pivoted = T, #output to be pivoted (true) or not (false)
                             unpivoted_with_bands = T,#output pivoted with bands (T) or with wavelength (false)
                             bands_HYP = read.csv("//10.0.1.243/projects/2023_HyperNet-POP/6-script/bands_HYP.csv")){
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
    hyp_wvl  <-ncdf4::ncvar_get(ncin, "wavelength")
    hyp_vaa <-ncdf4::ncvar_get(ncin, "viewing_azimuth_angle")
    hyp_paa   <-ncdf4::ncvar_get(ncin, "pointing_azimuth_angle")
    hyp_vza  <-ncdf4::ncvar_get(ncin, "viewing_zenith_angle")
    hyp_sza <- ncdf4::ncvar_get(ncin, "solar_zenith_angle")
    hyp_saa <- ncdf4::ncvar_get(ncin, "solar_azimuth_angle")
    
    complete_DB <- data.table::data.table(date = date,
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
    
    
    return(complete_DB)
    
  }) %>% bind_rows() %>%
    filter(B1000 < 6000 ) # max value @ 836 nm
  
  
  res2_long <- tidyr::pivot_longer(res, cols = starts_with("B"),
                                   names_to = "bands",
                                   values_to = "Ref",
                                   values_drop_na = TRUE) 
  res2_long <-left_join (res2_long,bands_HYP)
  res <- res2_long 
  
  #remove PAA
  
  res <- res  %>% select(-PAA) %>%
    mutate(time = substring(time,2,8)) %>%
    mutate(time = sub("\\s+$", "", gsub('(.{2})', '\\1:', time))) %>%
    mutate(time = substring(time,1,8))
    # mutate(SAA = correct_azimuth_angle(SAA)) %>%
    # mutate(VAA = correct_azimuth_angle(VAA)) %>%
    # mutate(SZA = abs(SZA)) %>%
    # mutate(VZA = abs(VZA))
    # 
  if(!output_pivoted){
    if(unpivoted_with_bands){
      not_pivoted <- res %>% select(-wvl) %>%
        tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
    }else{
      not_pivoted <- res %>% select(-bands) %>%
        tidyr::pivot_wider(names_from = "wvl", values_from = "Ref")
    }
    
    return(not_pivoted)
  }
  
  return(res)
}



#use res variable to make matchup with satellite image
#hypstar_angle_filter_and_mean ----
hypstar_angle_filter_and_mean <- function(dataset_res, #dataset pivoted of (date,time,VZA,VAA,PAA,B*) #res2_long oppure res
                            minVZA = min(unique(dataset_res$VZA)), #minimum View Zenith Angle
                            maxVZA = max(unique(dataset_res$VZA)), #maximum View Zenith Angle
                            minVAA = min(unique(dataset_res$VAA)), #minimum View Azimuth Angle
                            maxVAA = max(unique(dataset_res$VAA)), #maximum View Azimuth Angle
                            minSZA = min(unique(dataset_res$SZA)), #minimum Solar Zenith Angle
                            maxSZA = max(unique(dataset_res$SZA)), #maximum Solar Zenith Angle
                            minSAA = min(unique(dataset_res$SAA)), #minimum Solar Azimuth Angle
                            maxSAA = max(unique(dataset_res$SAA)), #maximum Solar Azimuth Angle
                            delVZA = c(""), #vector deleting View Zenith Angle inside dataset
                            delSZA = c(""), #vector deleting Solar Zenith Angle inside dataset
                            delVAA = c(""), #vector deleting View Azimuth Angle inside datasets
                            delSAA = c(""), #vector deleting Solar Azimuth Angle inside datasets
                            input_pivoted = T, #dataset_res is pivoted (true) or not (false) 
                            output_pivoted = T, #output to be pivoted (true) or not (false)
                            unpivoted_with_bands = T,#output pivoted with bands (T) or with wavelength (false)
                            bands_HYP = read.csv("//10.0.1.243/projects/2023_HyperNet-POP/6-script/bands_HYP.csv")
                            ){
  library(tidytable)

  if(!input_pivoted){
    res2_long <- tidyr::pivot_longer(dataset_res, cols = starts_with("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
    res2_long <-left_join (res2_long,bands_HYP)
    dataset_res <- res2_long  
  }
  
  filtered_res <- dataset_res %>% filter(VZA >= minVZA & VZA <= maxVZA) %>%
    filter(VAA >= minVAA & VAA <= maxVAA) %>%
    filter(SZA >= minSZA & SZA <= maxSZA) %>%
    filter(SAA >= minSAA & SAA <= maxSAA) %>%
    filter(!(VZA %in% delVZA)) %>%
    filter(!(VAA %in% delVAA)) %>%
    filter(!(SZA %in% delSZA)) %>%
    filter(!(SAA %in% delSAA)) %>%
    mutate(VZA_VAA_SZA_SAA = paste("(",VZA,":",VAA,":",SZA,":",SAA,")",sep="")) 
  
  VZA_VAA_SZA_SAA_used <- unique(filtered_res$VZA_VAA_SZA_SAA)
  
  matchup_selection_res <- filtered_res %>% 
    group_by(date,time,bands) %>%
    mutate(Ref_mean = mean(Ref)) %>%
    ungroup() %>%
    select(-Ref,-VZA,-VAA,-SZA,-SAA,-VZA_VAA_SZA_SAA) %>%
    rename(Ref = Ref_mean) %>%
    unique() %>%
    mutate(VZA_VAA_SZA_SAA_used = paste(VZA_VAA_SZA_SAA_used, collapse = ","))
  
  if(!output_pivoted){
    if(unpivoted_with_bands){
      not_pivoted <- matchup_selection_res %>% select(-wvl) %>%
        tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
    }else{
      not_pivoted <- matchup_selection_res %>% select(-bands) %>%
        tidyr::pivot_wider(names_from = "wvl", values_from = "Ref")
    }
    
    return(not_pivoted)
  }
  
  return(matchup_selection_res)
  
}

#hypstar_spectral_resampling_to_S2 ----
hypstar_spectral_resampling_to_S2 <- function(dataset_res,
                                        S2_sensor,#can be S2A,S2B,S2C
                                        input_pivoted = T, #dataset_res is pivoted (true) or not (false) 
                                        output_pivoted = T, #output to be pivoted (true) or not (false)
                                        unpivoted_with_bands = T,#output pivoted with bands (T) or with wavelength (false)
                                        bands_HYP = read.csv("//10.0.1.243/projects/2023_HyperNet-POP/6-script/bands_HYP.csv"),
                                        S2_bands = c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"),
                                        bands_S2 = read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/s2_bands.csv")){
  library(tidytable)
  
  if(!input_pivoted){
    res2_long <- tidyr::pivot_longer(dataset_res, cols = starts_with("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
    res2_long <-left_join (res2_long,bands_HYP)
    dataset_res <- res2_long  
  }
  
  source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/F01_Spectral_Resampling_V02.R")
  
  origin_band_centers <- sort(unique(bands_HYP$wvl))
  
  #Crea la matrice di ricampionamento
  SFR_S2 <- suppressWarnings(SRF_adapter(origin_band_centers = origin_band_centers,
                        mode = "SRF",
                        sensor = as.character(S2_sensor),
                        ref_folder = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/config_data")[,S2_bands])
  
  #Eseguo il ricampionamento
  unpivoted_res <- dataset_res %>% select(-wvl) %>%
    tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
  
  band_variables <- names(unpivoted_res)[grepl(x = names(unpivoted_res), pattern = "^B.*")]
  
  file_resampled <- dataframe_resampler(total_df = unpivoted_res, 
                                               band_variables = band_variables, #qui serve il nome che le bande hanno nel file originale
                                               SRF = SFR_S2) 
  
  res2_long <- tidyr::pivot_longer(file_resampled, cols = starts_with("B"),
                                   names_to = "bands",
                                   values_to = "Ref",
                                   values_drop_na = TRUE) 
  res2_long <-left_join (res2_long,bands_S2)
  file_resampled <- res2_long
  
  if(!output_pivoted){
    if(unpivoted_with_bands){
      not_pivoted <- file_resampled %>% select(-wvl) %>%
        tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
    }else{
      not_pivoted <- file_resampled %>% select(-bands) %>%
        tidyr::pivot_wider(names_from = "wvl", values_from = "Ref")
    }
    
    return(not_pivoted)
  }
  
  return(file_resampled)
  
}

#hypstar_resampling_to_PRISMA ----
hypstar_spectral_resampling_to_PRISMA <- function(dataset_res,
                                                  input_pivoted = T, #dataset_res is pivoted (true) or not (false)
                                                  output_pivoted = T, #output to be pivoted (true) or not (false)
                                                  unpivoted_with_bands = T,#output pivoted with bands (T) or with wavelength (false)
                                                  bands_HYP = read.csv("//10.0.1.243/projects/2023_HyperNet-POP/6-script/bands_HYP.csv"),
                                                  PRISMA_bands = paste0("B",1:151),#select which PRISMA_bands to choose (only hypstar bands for PRISMA)
                                                  bands_PRISMA = read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/PRISMA_bands.csv")){
  library(tidytable)
  
  if(!input_pivoted){
    res2_long <- tidyr::pivot_longer(dataset_res, cols = starts_with("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
    res2_long <-left_join (res2_long,bands_HYP)
    dataset_res <- res2_long  
  }
  
  source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/F01_Spectral_Resampling_V02.R")
  
  origin_band_centers <- sort(unique(bands_HYP$wvl))
  
  #Crea la matrice di ricampionamento
  SFR_PRISMA <- suppressWarnings(SRF_adapter(origin_band_centers = origin_band_centers,
                            mode = "FWHM",
                            sensor = "PRISMA",
                            ref_folder = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/config_data"))
  
  #Eseguo il ricampionamento
  unpivoted_res <- dataset_res %>% select(-wvl) %>%
    tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
  
  band_variables <- names(unpivoted_res)[grepl(x = names(unpivoted_res), pattern = "^B.*")]
  
  file_resampled <- dataframe_resampler(total_df = unpivoted_res, 
                                        band_variables = band_variables, #qui serve il nome che le bande hanno nel file originale
                                        SRF = SFR_PRISMA) 
  
  res2_long <- tidyr::pivot_longer(file_resampled, cols = starts_with("B"),
                                   names_to = "bands",
                                   values_to = "Ref",
                                   values_drop_na = TRUE) 
  res2_long <-left_join (res2_long,bands_PRISMA)
  file_resampled <- res2_long %>% na.omit()
  
  if(!output_pivoted){
    if(unpivoted_with_bands){
      not_pivoted <- file_resampled %>% select(-wvl) %>%
        tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
    }else{
      not_pivoted <- file_resampled %>% select(-bands) %>%
        tidyr::pivot_wider(names_from = "wvl", values_from = "Ref")
    }
    
    return(not_pivoted)
  }
  
  return(file_resampled)
  
}

#get_s2_data ----
get_s2_data <- function(date_start, 
                        date_end,
                        hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                        repo_root_path = "//10.0.1.243/Repository/nr_products/sentinel2/l2/T32TQQ/", #masked S2 data
                        output_pivoted = T, #output to be pivoted (true) or not (false)
                        unpivoted_with_bands = T,#output pivoted with bands (T) or with wavelength (false)
                        bands_S2 = read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/s2_bands.csv")){
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
      output <- mean(four_closest_pixels$value, na.rm = T)
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
    
    statistics_dataset <- data.table::data.table()
    
    statistics_pivoted <- statistics_dataset %>% mutate(bands = extract_names) %>%
      mutate(Ref = as.numeric(bands_statistics)) %>%
      mutate(date = as.Date(substring(basename(product_path),8,15),"%Y%m%d")) %>%
      mutate(time = substring(basename(product_path),16,22)) %>% 
      mutate(time = substring(time,2,8)) %>%
      mutate(time = sub("\\s+$", "", gsub('(.{2})', '\\1:', time))) %>%
      mutate(time = substring(time,1,8)) %>%
      pivot_wider(names_from = bands, values_from = Ref) %>%
      mutate(sensor = sensor_ID) %>%
      data.table::as.data.table()
    
    #statistics_dataset$bands <- extract_names
    
    #statistics_dataset$Ref <- as.numeric(bands_statistics)
    
    #here we want to extract also SZA, SAA, VZA, VAA
    
    return(statistics_pivoted)
  }
  
  all_datasets <- lapply(product_paths,get_statistics_dataset)
  
  res <- dplyr::bind_rows(all_datasets) %>%
    na.omit() #%>%
    # rowwise() %>%
    # mutate(SAA = correct_azimuth_angle(SAA)) %>%
    # ungroup() %>%
    # rowwise() %>%
    # mutate(VAA = correct_azimuth_angle(VAA)) %>%
    # ungroup() %>%
    # rowwise() %>%
    # mutate(SZA = abs(SZA)) %>%
    # ungroup() %>%
    # rowwise() %>%
    # mutate(VZA = abs(VZA)) %>%
    # ungroup()
    # 
  res2_long <- tidyr::pivot_longer(res, cols = starts_with("B"),
                                   names_to = "bands",
                                   values_to = "Ref",
                                   values_drop_na = TRUE) 
  res2_long <-left_join (res2_long,bands_S2)
  
  res <- res2_long
  
  if(!output_pivoted){
    if(unpivoted_with_bands){
      not_pivoted <- res %>% select(-wvl) %>%
        tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
    }else{
      not_pivoted <- res %>% select(-bands) %>%
        tidyr::pivot_wider(names_from = "wvl", values_from = "Ref")
    }
    
    return(not_pivoted)
  }

  return(res)
}

#s2_angle_filter ----
s2_angle_filter <- function(dataset_res, #dataset pivoted of (date,time,VZA,VAA,PAA,B*) #res2_long oppure res
                               minVZA = min(unique(dataset_res$VZA)), #minimum View Zenith Angle
                               maxVZA = max(unique(dataset_res$VZA)), #maximum View Zenith Angle
                               minVAA = min(unique(dataset_res$VAA)), #minimum View Azimuth Angle
                               maxVAA = max(unique(dataset_res$VAA)), #maximum View Azimuth Angle
                               minSZA = min(unique(dataset_res$SZA)), #minimum Solar Zenith Angle
                               maxSZA = max(unique(dataset_res$SZA)), #maximum Solar Zenith Angle
                               minSAA = min(unique(dataset_res$SAA)), #minimum Solar Azimuth Angle
                               maxSAA = max(unique(dataset_res$SAA)), #maximum Solar Azimuth Angle
                               input_pivoted = T, #dataset_res is pivoted (true) or not (false) 
                               output_pivoted = T, #output to be pivoted (true) or not (false)
                               unpivoted_with_bands = T,#output pivoted with bands (T) or with wavelength (false)
                               bands_S2 = read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/s2_bands.csv") #path bands_S2
                               ){
  if(!input_pivoted){
    res2_long <- tidyr::pivot_longer(dataset_res, cols = starts_with("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
    res2_long <-left_join(res2_long,bands_S2)
    dataset_res <- res2_long  
  }
  
  filtered_res <- dataset_res %>% filter(VZA >= minVZA & VZA <= maxVZA) %>%
    filter(VAA >= minVAA & VAA <= maxVAA) %>%
    filter(SZA >= minSZA & SZA <= maxSZA) %>%
    filter(SAA >= minSAA & SAA <= maxSAA) #%>%
    #filter(!(VZA %in% delVZA)) %>%
    #filter(!(VAA %in% delVAA)) %>%
    #filter(!(SZA %in% delSZA)) %>%
    #filter(!(SAA %in% delSAA)) %>%
    #mutate(VZA_VAA_SZA_SAA = paste("(",VZA,":",VAA,":",SZA,":",SAA,")",sep="")) 
  if(!output_pivoted){
    if(unpivoted_with_bands){
      not_pivoted <- filtered_res %>% select(-wvl) %>%
        tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
    }else{
      not_pivoted <- filtered_res %>% select(-bands) %>%
        tidyr::pivot_wider(names_from = "wvl", values_from = "Ref")
    }
    
    return(not_pivoted)
  }
  
  return(filtered_res)
}

#get_PRISMA_data ----
get_PRISMA_data <- function(date_start, 
                            date_end,
                            PRISMA_level_product = "L1",
                            processing_level = "Atcor_regrid_crop_smooth",
                            hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                            repo_root_path = "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/", #PRISMA data
                            output_pivoted = T, #output to be pivoted (true) or not (false)
                            unpivoted_with_bands = T,#output pivoted with bands (T) or with wavelength (false)
                            bands_PRISMA = read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/PRISMA_bands.csv")){
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
      output <- mean(four_closest_pixels$value, na.rm = T)
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
    
    bands_statistics <- lapply(bands, get_four_closest_pixels_statistic)
    
    statistics_dataset <- data.table::data.table()
    
    if(grepl("OFFL",basename(product_path))){
      OFFL_offset <- 5
    }else{
      OFFL_offset <- 0
    }
    if(grepl("STD",basename(product_path))){
      STD_offset <- 4
    }else{
      STD_offset <- 0
    }
    
    read_angle_and_time <- read.csv(paste0(product_path,"/ATCOR/all_angles_file.csv"))
    
    statistics_pivoted <- statistics_dataset %>% mutate(bands = paste0("B",1:230)) %>%
      mutate(Ref = as.numeric(bands_statistics)) %>%
      mutate(date = as.Date(substring(read_angle_and_time$date, 1,10),"%Y-%m-%d")) %>%
      mutate(time = substring(read_angle_and_time$date, 12,19)) %>% 
      pivot_wider(names_from = bands, values_from = Ref) %>%
      mutate(sensor = sensor_ID) %>%
      mutate(SAA = read_angle_and_time$sunaz) %>%
      mutate(SZA = read_angle_and_time$sunzen) %>%
      mutate(VZA = read_angle_and_time$sensor_zenith) %>%
      mutate(VAA = read_angle_and_time$sensor_azimuth) %>%
      data.table::as.data.table()
    
    #statistics_dataset$bands <- extract_names
    
    #statistics_dataset$Ref <- as.numeric(bands_statistics)
    
    #here we want to extract also SZA, SAA, VZA, VAA
    
    return(statistics_pivoted)
  }
  
  all_datasets <- lapply(product_paths,get_statistics_dataset)
  
  res <- dplyr::bind_rows(all_datasets) %>%
    na.omit() #%>%
    # rowwise() %>%
    # mutate(SAA = correct_azimuth_angle(SAA)) %>%
    # ungroup() %>%
    # rowwise() %>%
    # mutate(VAA = correct_azimuth_angle(VAA)) %>%
    # ungroup() %>%
    # rowwise() %>%
    # mutate(SZA = abs(SZA)) %>%
    # ungroup() %>%
    # rowwise() %>%
    # mutate(VZA = abs(VZA)) %>%
    # ungroup()
    # 
  res2_long <- tidyr::pivot_longer(res, cols = starts_with("B"),
                                   names_to = "bands",
                                   values_to = "Ref",
                                   values_drop_na = TRUE) 
  res2_long <-left_join (res2_long,bands_PRISMA)
  
  res <- res2_long
  
  if(!output_pivoted){
    if(unpivoted_with_bands){
      not_pivoted <- res %>% select(-wvl) %>%
        tidyr::pivot_wider(names_from = "bands", values_from = "Ref")
    }else{
      not_pivoted <- res %>% select(-bands) %>%
        tidyr::pivot_wider(names_from = "wvl", values_from = "Ref")
    }
    
    return(not_pivoted)
  }
  
  return(res)
  
  
}

#matchup_function_to_S2 ----
matchup_function_to_S2 <- function(s2_data, 
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
                             
                             input_pivoted = T, #dataset_res is pivoted (true) or not (false) 
                             output_pivoted = T, #output to be pivoted (true) or not (false)
                             unpivoted_with_bands = T,#output pivoted with bands (T) or with wavelength (false)
                             bands_HYP = read.csv("//10.0.1.243/projects/2023_HyperNet-POP/6-script/bands_HYP.csv"),
                             bands_S2 = read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/s2_bands.csv")
                             ){
  if(!input_pivoted){
    res2_long <- tidyr::pivot_longer(s2_data, cols = starts_with("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
    res2_long <-left_join(res2_long,bands_S2)
    s2_data <- res2_long  
    
    res2_long <- tidyr::pivot_longer(hypstar_data, cols = starts_with("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
    res2_long <-left_join(res2_long,bands_HYP)
    hypstar_data <- res2_long 
  }
  
  #to remove specific sensor view angles
  hypstar_data <- hypstar_data %>%     
    filter(!(VZA %in% delVZA)) %>%
    filter(!(VAA %in% delVAA))
    
  #check spectral resampling of hypstar
  hypstar_bands_number <- length(unique(hypstar_data$bands))
  if(hypstar_bands_number > 13){
    print("mancato ricampionamento, viene eseguito")
    hypstar_resampled_to_S2A <- hypstar_spectral_resampling_to_S2(hypstar_data,
                                                           S2_sensor = "S2A", 
                                                           input_pivoted = T) %>%
      mutate(sensor = "S2A")
    
    hypstar_resampled_to_S2B <- hypstar_spectral_resampling_to_S2(hypstar_data,
                                                                  S2_sensor = "S2B", 
                                                                  input_pivoted = T) %>%
      mutate(sensor = "S2B")
    
    hypstar_resampled_to_S2C <- hypstar_spectral_resampling_to_S2(hypstar_data,
                                                                  S2_sensor = "S2C", 
                                                                  input_pivoted = T) %>%
      mutate(sensor = "S2C")
    
    hypstar_resampled <- hypstar_resampled_to_S2A %>% bind_rows(hypstar_resampled_to_S2B) %>%
      bind_rows(hypstar_resampled_to_S2C)
    
  }else{
    hypstar_resampled <- hypstar_data
  }
  
  #change col names not to confuse S2 B# and hypstar B#
  s2_col_renamed <- s2_data %>% rename_with(~paste0(.,"_s2"),-c(bands,date,sensor,wvl)) 
  
  hypstar_col_renamed <- hypstar_resampled %>% rename_with(~paste0(.,"_hyp"),-c(bands,date,sensor,wvl)) 
  
    
  library(tidytable)
  
  
  #DATE PROBLEM
  if(only_same_date){
    #only_same_date
    date_output <- s2_col_renamed %>% left_join(hypstar_col_renamed, by = c("bands","date","sensor","wvl")) %>%
      na.omit()
     
  }else{
    #closest date 
    print("closest date will have the S2 date as 'date' variable")
    s2_with_date <- s2_col_renamed %>% rename(date_s2 = date)
    hypstar_with_date <- hypstar_col_renamed %>% rename(date_hyp = date)
    date_output <- s2_with_date %>% left_join(hypstar_with_date, by = c("bands","sensor","wvl")) %>%
      mutate(date_diff = abs(difftime(date_s2, date_hyp))) %>%
      group_by(sensor,bands) %>%
      mutate(lowest_date_diff = min(date_diff)) %>%
      ungroup() %>%
      filter(date_diff == lowest_date_diff) %>%
      select(-date_diff,-lowest_date_diff) %>%
      select(-date_hyp) %>%
      rename(date = date_s2) %>%
      na.omit()
  }
  
  #SUN TIME OR SUN ANGLE PROBLEM
  if(time_or_sun_angle){
    
    #create time variables and time difference
    time_output <- date_output %>%
      mutate(date_time_s2 = as.POSIXct(paste(date,time_s2))) %>%
      mutate(date_time_hyp = as.POSIXct(paste(date,time_hyp))) %>%
      mutate(time_diff = abs(difftime(date_time_s2, date_time_hyp, units = "sec")))%>%
      select(-date_time_s2,-date_time_hyp) 
    
    if(just_closest_time_or_sun_angle){
      #only closest time
      time_or_sun_angle_output <- time_output %>% 
        group_by(date,bands,sensor,wvl) %>%
        mutate(lowest_time_diff = min(time_diff)) %>%
        ungroup() %>%
        filter(time_diff == lowest_time_diff) %>%
        select(-time_diff,-lowest_time_diff) 
    }else{
      #take a window of times and mean
      
      #here we need to use the geosphere package to compute the centroid of a spherical polygon
      convert_angles_for_geosphere <- time_output %>%
        mutate(SAA_hyp = ifelse(SAA_hyp > 180, 180-SAA_hyp, SAA_hyp )) 
      
      reconvert_for_us <- convert_angles_for_geosphere %>%
        filter(time_diff < range_time_or_sun_angle/2*3600) %>%
        select(-time_diff) %>%
        select(-time_hyp) %>%
        group_by(date,bands,sensor,wvl,VZA_hyp,VAA_hyp) %>%
        mutate(count = row_number()) %>%
        mutate(mean_Ref_hyp = mean(Ref_hyp)) %>%
        mutate(mean_SZA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[2]) %>%
        mutate(mean_SAA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[1]) %>%
        ungroup() %>%
        select(-count) %>%
        select(-Ref_hyp) %>%
        select(-SZA_hyp) %>%
        select(-SAA_hyp) %>%
        #select(-time_diff)%>%
        unique() 
      
      time_or_sun_angle_output <- reconvert_for_us %>% 
        mutate(mean_SAA_hyp = ifelse(mean_SAA_hyp < 0, 180 - mean_SAA_hyp, mean_SAA_hyp))
      
      # if(weighted_mean){
      #   #create a weighted mean over times
      #   time_or_sun_output <- time_output %>%
      #     filter(time_diff < range_time_or_sun_angle/2*3600) %>%
      #     #select(-time_diff) %>%
      #     select(-time_hyp) %>%
      #     mutate(num_Ref = as.numeric(time_diff)*Ref_hyp) %>%
      #     mutate(num_SZA = as.numeric(time_diff)*SZA_hyp) %>%
      #     mutate(num_SAA = as.numeric(time_diff)*SAA_hyp) %>%
      #     group_by(date,bands,sensor,VZA_hyp,VAA_hyp) %>%
      #     mutate(count = row_number()) %>%
      #     mutate(total_den = sum(as.numeric(time_diff))) %>%
      #     mutate(total_num_Ref = sum(num_Ref)) %>%
      #     mutate(total_num_SZA = sum(num_SZA)) %>%
      #     mutate(total_num_SAA = sum(num_SAA)) %>%
      #     ungroup() %>%
      #     mutate(mean_Ref_hyp = total_num_Ref/total_den) %>%
      #     mutate(mean_SZA_hyp = total_num_SZA/total_den) %>%
      #     mutate(mean_SAA_hyp = total_num_SAA/total_den) %>%
      #     select(-count) %>%
      #     select(-Ref_hyp,-SZA_hyp,-SAA_hyp) %>%
      #     select(-num_Ref,-num_SZA,-num_SAA) %>%
      #     select(-total_den,-total_num_Ref,-total_num_SZA,-total_num_SAA) %>%
      #     select(-time_diff) %>%
      #     unique()
      # }else{#execute mean over all times in window without weight
      #   
      # }
      
    }
  }else{
    #create time variables and angle difference
    sun_angle_output <- date_output %>%
      mutate(time_s2 = substring(time_s2,2,8)) %>%
      mutate(time_s2 = sub("\\s+$", "", gsub('(.{2})', '\\1:', time_s2))) %>%
      mutate(time_s2 = substring(time_s2,1,8)) %>%
      mutate(time_hyp = substring(time_hyp,2,8)) %>%
      mutate(time_hyp = sub("\\s+$", "", gsub('(.{2})', '\\1:', time_hyp))) %>%
      mutate(time_hyp = substring(time_hyp,1,8)) %>%
      mutate(angle_diff = 180/pi*acos( sin(SAA_s2*pi/180)*sin(SAA_hyp*pi/180) + cos(SAA_s2*pi/180) * cos(SAA_hyp*pi/180) * cos( SZA_s2*pi/180-SZA_hyp*pi/180 ))) 
    #https://en.wikipedia.org/wiki/Great-circle_distance#cite_ref-2
      
    if(just_closest_time_or_sun_angle){
      #only closest angle
      time_or_sun_angle_output <- sun_angle_output %>% 
        group_by(date,bands,wvl,sensor) %>%
        mutate(lowest_angle_diff = min(angle_diff)) %>%
        ungroup() %>%
        filter(angle_diff == lowest_angle_diff) %>%
        select(-angle_diff,-lowest_angle_diff) 
    }else{
      #take a window of angles and mean
      
      convert_angles_for_geosphere <- sun_angle_output %>%
        mutate(SAA_hyp = ifelse(SAA_hyp > 180, 180-SAA_hyp, SAA_hyp )) 
      
      reconvert_for_us <- convert_angles_for_geosphere %>%
        filter(angle_diff < range_time_or_sun_angle/2) %>%
        select(-angle_diff) %>%
        select(-time_hyp) %>%
        group_by(date,bands,sensor,wvl,VZA_hyp,VAA_hyp) %>%
        mutate(count = row_number()) %>%
        mutate(mean_Ref_hyp = mean(Ref_hyp)) %>%
        mutate(mean_SZA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[2]) %>%
        mutate(mean_SAA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[1]) %>%
        ungroup() %>%
        select(-count) %>%
        select(-Ref_hyp) %>%
        select(-SZA_hyp) %>%
        select(-SAA_hyp) %>%
        #select(-angle_diff)%>%
        unique() 
      
      time_or_sun_angle_output <- reconvert_for_us %>% 
        mutate(mean_SAA_hyp = ifelse(mean_SAA_hyp < 0, 180 - mean_SAA_hyp, mean_SAA_hyp))
      
    }
  }
  
  #VIEW ANGLE PROBLEM
  med_view_angle_output <- time_or_sun_angle_output %>%
    mutate(angle_diff = 180/pi*acos( sin(VAA_s2*pi/180)*sin(VAA_hyp*pi/180) + cos(VAA_s2*pi/180) * cos(VAA_hyp*pi/180) * cos( VZA_s2*pi/180-VZA_hyp*pi/180 ))) 
  #https://en.wikipedia.org/wiki/Great-circle_distance#cite_ref-2
  
  if(just_closest_view_angle){
    #only closest angle
    view_angle_output <- med_view_angle_output %>% 
      group_by(date,bands,wvl,sensor) %>%
      mutate(lowest_angle_diff = min(angle_diff)) %>%
      ungroup() %>%
      filter(angle_diff == lowest_angle_diff) %>%
      select(-angle_diff,-lowest_angle_diff) 
  }else{
    #take a window of angles and mean
    
    convert_angles_for_geosphere <- med_view_angle_output %>%
      mutate(VAA_hyp = ifelse(VAA_hyp > 180, VAA_hyp - 360, VAA_hyp )) 
    
    reconvert_for_us <- convert_angles_for_geosphere %>%
      filter(angle_diff < range_view_angle/2) %>%
      select(-angle_diff) %>%
      select(-time_hyp) %>%
      group_by(date,bands,wvl,sensor) %>%
      mutate(count = row_number()) %>%
      mutate(mean_Ref_hyp = mean(Ref_hyp)) %>%
      mutate(mean_VZA_hyp = geosphere::centroid(cbind(VAA_hyp,VZA_hyp))[2]) %>%
      mutate(mean_VAA_hyp = geosphere::centroid(cbind(VAA_hyp,VZA_hyp))[1]) %>%
      ungroup() %>%
      select(-count) %>%
      select(-Ref_hyp) %>%
      select(-VZA_hyp) %>%
      select(-VAA_hyp) %>%
      #select(-angle_diff)%>%
      unique() 
    
    view_angle_output <- reconvert_for_us %>% 
      mutate(mean_VAA_hyp = ifelse(mean_VAA_hyp < 0, 360 + mean_VAA_hyp, mean_VAA_hyp))
    
  }
  
  return(view_angle_output)
}

#matchup_function_to_PRISMA ----
matchup_function_to_PRISMA <- function(PRISMA_data, 
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
                                       
                                       input_pivoted = T, #dataset_res is pivoted (true) or not (false) 
                                       output_pivoted = T, #output to be pivoted (true) or not (false)
                                       unpivoted_with_bands = T,#output pivoted with bands (T) or with wavelength (false)
                                       bands_HYP = read.csv("//10.0.1.243/projects/2023_HyperNet-POP/6-script/bands_HYP.csv"),
                                       bands_PRISMA = read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/PRISMA_bands.csv")){
  
  if(!input_pivoted){
    res2_long <- tidyr::pivot_longer(PRISMA_data, cols = starts_with("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
    res2_long <-left_join(res2_long,bands_PRISMA)
    PRISMA_data <- res2_long  
    
    res2_long <- tidyr::pivot_longer(hypstar_data, cols = starts_with("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
    res2_long <-left_join(res2_long,bands_HYP)
    hypstar_data <- res2_long 
  }
  
  #to remove specific sensor view angles
  hypstar_data <- hypstar_data %>%     
    filter(!(VZA %in% delVZA)) %>%
    filter(!(VAA %in% delVAA))
  
  #check spectral resampling of hypstar
  hypstar_bands_number <- length(unique(hypstar_data$bands))
  if(hypstar_bands_number > 151){
    print("mancato ricampionamento, viene eseguito")
    hypstar_resampled <- hypstar_spectral_resampling_to_PRISMA(hypstar_data,
                                                               input_pivoted = T) %>%
      mutate(sensor = "PRISMA")
    
  }else{
    hypstar_resampled <- hypstar_data
  }
  
  #change col names not to confuse PRISMA B# and hypstar B#
  PRISMA_col_renamed <- PRISMA_data %>% rename_with(~paste0(.,"_prs"),-c(bands,date,sensor,wvl)) 
  
  hypstar_col_renamed <- hypstar_resampled %>% rename_with(~paste0(.,"_hyp"),-c(bands,date,sensor,wvl)) 
  
  
  library(tidytable)
  
  
  #DATE PROBLEM
  if(only_same_date){
    #only_same_date
    date_output <- PRISMA_col_renamed %>% left_join(hypstar_col_renamed, by = c("bands","date","sensor","wvl")) %>%
      na.omit()
    
  }else{
    #closest date 
    print("closest date will have the S2 date as 'date' variable")
    PRISMA_with_date <- PRISMA_col_renamed %>% rename(date_prs = date)
    hypstar_with_date <- hypstar_col_renamed %>% rename(date_hyp = date)
    date_output <- PRISMA_with_date %>% left_join(hypstar_with_date, by = c("bands","sensor","wvl")) %>%
      mutate(date_diff = abs(difftime(date_prs, date_hyp))) %>%
      group_by(sensor,bands) %>%
      mutate(lowest_date_diff = min(date_diff)) %>%
      ungroup() %>%
      filter(date_diff == lowest_date_diff) %>%
      select(-date_diff,-lowest_date_diff) %>%
      select(-date_hyp) %>%
      rename(date = date_prs) %>%
      na.omit()
  }
  
  #SUN TIME OR SUN ANGLE PROBLEM
  if(time_or_sun_angle){
    
    #create time variables and time difference
    time_output <- date_output %>%
      mutate(date_time_prs = as.POSIXct(paste(date,time_prs))) %>%
      mutate(date_time_hyp = as.POSIXct(paste(date,time_hyp))) %>%
      mutate(time_diff = abs(difftime(date_time_prs, date_time_hyp, units = "sec")))%>%
      select(-date_time_prs,-date_time_hyp) 
    
    if(just_closest_time_or_sun_angle){
      #only closest time
      time_or_sun_angle_output <- time_output %>% 
        group_by(date,bands,sensor,wvl) %>%
        mutate(lowest_time_diff = min(time_diff)) %>%
        ungroup() %>%
        filter(time_diff == lowest_time_diff) %>%
        select(-time_diff,-lowest_time_diff) 
    }else{
      #take a window of times and mean
      
      #here we need to use the geosphere package to compute the centroid of a spherical polygon
      convert_angles_for_geosphere <- time_output %>%
        mutate(SAA_hyp = ifelse(SAA_hyp > 180, SAA_hyp-360, SAA_hyp )) 
      
      reconvert_for_us <- convert_angles_for_geosphere %>%
        filter(time_diff < range_time_or_sun_angle/2*3600) %>%
        select(-time_diff) %>%
        select(-time_hyp) %>%
        group_by(date,bands,sensor,wvl,VZA_hyp,VAA_hyp) %>%
        mutate(count = row_number()) %>%
        mutate(mean_Ref_hyp = mean(Ref_hyp)) %>%
        mutate(mean_SZA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[2]) %>%
        mutate(mean_SAA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[1]) %>%
        ungroup() %>%
        select(-count) %>%
        select(-Ref_hyp) %>%
        select(-SZA_hyp) %>%
        select(-SAA_hyp) %>%
        #select(-time_diff)%>%
        unique() 
      
      time_or_sun_angle_output <- reconvert_for_us %>% 
        mutate(mean_SAA_hyp = ifelse(mean_SAA_hyp < 0, 360 + mean_SAA_hyp, mean_SAA_hyp))
      
      # if(weighted_mean){
      #   #create a weighted mean over times
      #   time_or_sun_output <- time_output %>%
      #     filter(time_diff < range_time_or_sun_angle/2*3600) %>%
      #     #select(-time_diff) %>%
      #     select(-time_hyp) %>%
      #     mutate(num_Ref = as.numeric(time_diff)*Ref_hyp) %>%
      #     mutate(num_SZA = as.numeric(time_diff)*SZA_hyp) %>%
      #     mutate(num_SAA = as.numeric(time_diff)*SAA_hyp) %>%
      #     group_by(date,bands,sensor,VZA_hyp,VAA_hyp) %>%
      #     mutate(count = row_number()) %>%
      #     mutate(total_den = sum(as.numeric(time_diff))) %>%
      #     mutate(total_num_Ref = sum(num_Ref)) %>%
      #     mutate(total_num_SZA = sum(num_SZA)) %>%
      #     mutate(total_num_SAA = sum(num_SAA)) %>%
      #     ungroup() %>%
      #     mutate(mean_Ref_hyp = total_num_Ref/total_den) %>%
      #     mutate(mean_SZA_hyp = total_num_SZA/total_den) %>%
      #     mutate(mean_SAA_hyp = total_num_SAA/total_den) %>%
      #     select(-count) %>%
      #     select(-Ref_hyp,-SZA_hyp,-SAA_hyp) %>%
      #     select(-num_Ref,-num_SZA,-num_SAA) %>%
      #     select(-total_den,-total_num_Ref,-total_num_SZA,-total_num_SAA) %>%
      #     select(-time_diff) %>%
      #     unique()
      # }else{#execute mean over all times in window without weight
      #   
      # }
      
    }
  }else{
    #create time variables and angle difference
    sun_angle_output <- date_output %>%
      mutate(time_prs = substring(time_prs,2,8)) %>%
      mutate(time_prs = sub("\\s+$", "", gsub('(.{2})', '\\1:', time_prs))) %>%
      mutate(time_prs = substring(time_prs,1,8)) %>%
      mutate(time_hyp = substring(time_hyp,2,8)) %>%
      mutate(time_hyp = sub("\\s+$", "", gsub('(.{2})', '\\1:', time_hyp))) %>%
      mutate(time_hyp = substring(time_hyp,1,8)) %>%
      mutate(angle_diff = 180/pi*acos( sin(SAA_prs*pi/180)*sin(SAA_hyp*pi/180) + cos(SAA_prs*pi/180) * cos(SAA_hyp*pi/180) * cos( SZA_prs*pi/180-SZA_hyp*pi/180 ))) 
    #https://en.wikipedia.org/wiki/Great-circle_distance#cite_ref-2
    
    if(just_closest_time_or_sun_angle){
      #only closest angle
      time_or_sun_angle_output <- sun_angle_output %>% 
        group_by(date,bands,wvl,sensor) %>%
        mutate(lowest_angle_diff = min(angle_diff)) %>%
        ungroup() %>%
        filter(angle_diff == lowest_angle_diff) %>%
        select(-angle_diff,-lowest_angle_diff) 
    }else{
      #take a window of angles and mean
      
      convert_angles_for_geosphere <- sun_angle_output %>%
        mutate(SAA_hyp = ifelse(SAA_hyp > 180, SAA_hyp-360, SAA_hyp )) 
      
      reconvert_for_us <- convert_angles_for_geosphere %>%
        filter(angle_diff < range_time_or_sun_angle/2) %>%
        select(-angle_diff) %>%
        select(-time_hyp) %>%
        group_by(date,bands,sensor,wvl,VZA_hyp,VAA_hyp) %>%
        mutate(count = row_number()) %>%
        mutate(mean_Ref_hyp = mean(Ref_hyp)) %>%
        mutate(mean_SZA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[2]) %>%
        mutate(mean_SAA_hyp = geosphere::centroid(cbind(SAA_hyp,SZA_hyp))[1]) %>%
        ungroup() %>%
        select(-count) %>%
        select(-Ref_hyp) %>%
        select(-SZA_hyp) %>%
        select(-SAA_hyp) %>%
        #select(-angle_diff)%>%
        unique() 
      
      time_or_sun_angle_output <- reconvert_for_us %>% 
        mutate(mean_SAA_hyp = ifelse(mean_SAA_hyp < 0, 360 + mean_SAA_hyp, mean_SAA_hyp))
      
    }
  }
  
  #VIEW ANGLE PROBLEM
  med_view_angle_output <- time_or_sun_angle_output %>%
    mutate(angle_diff = 180/pi*acos( sin(VAA_s2*pi/180)*sin(VAA_hyp*pi/180) + cos(VAA_s2*pi/180) * cos(VAA_hyp*pi/180) * cos( VZA_s2*pi/180-VZA_hyp*pi/180 ))) 
  #https://en.wikipedia.org/wiki/Great-circle_distance#cite_ref-2
  
  if(just_closest_view_angle){
    #only closest angle
    view_angle_output <- med_view_angle_output %>% 
      group_by(date,bands,wvl,sensor) %>%
      mutate(lowest_angle_diff = min(angle_diff)) %>%
      ungroup() %>%
      filter(angle_diff == lowest_angle_diff) %>%
      select(-angle_diff,-lowest_angle_diff) 
  }else{
    #take a window of angles and mean
    
    convert_angles_for_geosphere <- med_view_angle_output %>%
      mutate(VAA_hyp = ifelse(VAA_hyp > 180, VAA_hyp-360, VAA_hyp )) 
    
    reconvert_for_us <- convert_angles_for_geosphere %>%
      filter(angle_diff < range_view_angle/2) %>%
      select(-angle_diff) %>%
      select(-time_hyp) %>%
      group_by(date,bands,wvl,sensor) %>%
      mutate(count = row_number()) %>%
      mutate(mean_Ref_hyp = mean(Ref_hyp)) %>%
      mutate(mean_VZA_hyp = geosphere::centroid(cbind(VAA_hyp,VZA_hyp))[2]) %>%
      mutate(mean_VAA_hyp = geosphere::centroid(cbind(VAA_hyp,VZA_hyp))[1]) %>%
      ungroup() %>%
      select(-count) %>%
      select(-Ref_hyp) %>%
      select(-VZA_hyp) %>%
      select(-VAA_hyp) %>%
      #select(-angle_diff)%>%
      unique() 
    
    view_angle_output <- reconvert_for_us %>% 
      mutate(mean_VAA_hyp = ifelse(mean_VAA_hyp < 0, 360 + mean_VAA_hyp, mean_VAA_hyp))
    
  }
  
  return(view_angle_output)
  
}






correct_azimuth_angle <- function(angle){
  if(angle < 0){
    return(180-angle)
  }else{
    return(angle)
  }
}



















