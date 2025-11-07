#Create match-up

get_hypstar_data <- function(home_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/processed",
                             date_start,
                             date_end,
                             pivoted = T,
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
  
  if(pivoted){
    res2_long <- tidyr::pivot_longer(res, cols = starts_with("B"),
                                     names_to = "bands",
                                     values_to = "Ref",
                                     values_drop_na = TRUE) 
    res2_long <-left_join (res2_long,bands_HYP)
    res <- res2_long 
  }
  
  return(res)
}



#use res variable to make matchup with satellite image

hypstar_select_measures_and_mean <- function(dataset_res, #dataset pivoted of (date,time,VZA,VAA,PAA,B*) #res2_long oppure res
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
    mutate(VZA_VAA_SZA_SAA = paste("(",VZA,";",VAA,";",SZA,";",SAA,")",sep="")) 
  
  VZA_VAA_SZA_SAA_used <- unique(filtered_res$VZA_VAA_SZA_SAA)
  
  matchup_selection_res <- filtered_res %>% 
    group_by(date,time,bands) %>%
    mutate(Ref_mean = mean(Ref)) %>%
    ungroup() %>%
    select(-Ref,-VZA,-VAA,-PAA,-SZA,-SAA,-VZA_VAA_SZA_SAA) %>%
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

get_s2_data <- function(hypstar_point_path = "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/VECTOR/HYPSTAR_station/Posizione_punto/JSIT_Hypstar_EPSG32632.shp",
                        repo_root_path = "//10.0.1.243/Repository/nr_products/sentinel2/l2/T32TQQ/", #masked S2 data
                        date_start, 
                        date_end,
                        output_pivoted = T){
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
      data.table::as.data.table()
    
    #statistics_dataset$bands <- extract_names
    
    #statistics_dataset$Ref <- as.numeric(bands_statistics)
    
    #here we want to extract also SZA, SAA, VZA, VAA
    
    return(statistics_pivoted)
  }
  
  all_datasets <- lapply(product_paths,get_statistics_dataset)
  
  output <- dplyr::bind_rows(all_datasets)   
  
  if(output_pivoted){
    output <- output %>% tidyr::pivot_longer(cols = starts_with("B"),
                                             names_to = "bands",
                                             values_to = "Ref") 
  }

  return(output)
}


#s2_select_measures_and_mean <- function()



