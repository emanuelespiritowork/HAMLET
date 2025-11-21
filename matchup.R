#Create match-up
#get_hypstar_data ----
get_hypstar_data <- function(home_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/processed",
                             date_start,
                             date_end,
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
  
  res <- res  %>% select(-PAA)
  
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

#get_s2_data ----
get_s2_data <- function(hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                        repo_root_path = "//10.0.1.243/Repository/nr_products/sentinel2/l2/T32TQQ/", #masked S2 data
                        date_start, 
                        date_end,
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
    na.omit()
  
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

#matchup_function ----
matchup_function_to_S2 <- function(s2_data, 
                             hypstar_data,
                             only_same_date = T,#only same date comparison accepted (true) or closest date (false)
                             just_closest_time = T,#closest time (true) or all times in a hour range (false)
                             hour_range = NULL,#if just_closest_time = F, here we need a range of hours. -hour_range/2,x,+hour_range/2
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
  s2_col_renamed <- s2_data %>% rename_with(~paste0(.,"_s2"),-c(bands,date,sensor)) %>%
    select(-wvl_s2) 
  
  hypstar_col_renamed <- hypstar_resampled %>% rename_with(~paste0(.,"_hyp"),-c(bands,date,sensor)) %>%
    select(-wvl_hyp)
  
    
  library(tidytable)
  
  
  #TIME PROBLEM
  if(only_same_date){
    merge_by_date_and_bands <- s2_col_renamed %>% left_join(hypstar_col_renamed, by = c("bands","date","sensor")) %>%
      mutate(time_s2 = substring(time_s2,2,8)) %>%
      mutate(time_s2 = sub("\\s+$", "", gsub('(.{2})', '\\1:', time_s2))) %>%
      mutate(time_s2 = substring(time_s2,1,8)) %>%
      mutate(time_hyp = substring(time_hyp,2,8)) %>%
      mutate(time_hyp = sub("\\s+$", "", gsub('(.{2})', '\\1:', time_hyp))) %>%
      mutate(time_hyp = substring(time_hyp,1,8)) %>%
      mutate(date_time_s2 = as.POSIXct(paste(date,time_s2))) %>%
      mutate(date_time_hyp = as.POSIXct(paste(date,time_hyp))) %>%
      mutate(time_diff = abs(difftime(date_time_s2, date_time_hyp, units = "sec")))%>%
      select(-date_time_s2,-date_time_hyp)
    
    if(just_closest_time){
      merge_closest_time <- merge_by_date_and_bands %>% 
        group_by(date,bands,sensor) %>%
        mutate(lowest_time_diff = min(time_diff)) %>%
        ungroup() %>%
        filter(time_diff == lowest_time_diff) %>%
        select(-time_diff,-lowest_time_diff) 
    }else{
      #execute mean over all times without weight
      merge_closest_time <- merge_by_date_and_bands %>%
        filter(time_diff < hour_range/2*3600) %>%
        select(-time_diff) %>%
        group_by(date,bands,sensor,VZA_hyp,VAA_hyp) %>%
        mutate(count = row_number()) %>%
        mutate(mean_Ref_hyp = mean(Ref_hyp)) %>%
        mutate(mean_SZA_hyp = mean(SZA_hyp)) %>%
        mutate(mean_SAA_hyp  = mean(SAA_hyp)) %>%
        ungroup() %>%
        select(-count) %>%
        select(-Ref_hyp) %>%
        select(-SZA_hyp) %>%
        select(-SAA_hyp) %>%
        select(-time_hyp) %>%
        #select(-time_diff)%>%
        unique() 
    }
     
  }
  
  #ANGLE PROBLEM
   
  
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
}




