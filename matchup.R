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

get_s2_data <- function(input_path = "//10.0.1.243/Repository/nr_products/sentinel2/l2/T32TQQ",
                        input_polygon,
                        output_file,
                        start_date,
                        end_date){
  terra::terraOptions(progress=0)
  
  # set arguments
  sp_res <- 20
  progress <- TRUE
  
  if(!(sp_res %in% c(10,20))){
    stop("Analysis spatial resolution should be one of the followings: '10', '20'.")
  }
  
  ptm <- proc.time()
  message(paste(c("Query started at: "), Sys.time(), sep=""))
  
  # list imagery in input folder
  f_list_full <- sort(list.files(path=input_path, pattern = "0m.vrt", recursive=TRUE, full.names = TRUE))
  
  # subset based on date
  f_list_full_date <- substr(basename(f_list_full), start = 8, stop = 15)
  start_fdate <- strftime(as.Date(start_date, origin="1970-01-01", tz="GMT"), format = "%Y%m%d", tz = "GMT", usetz = FALSE)
  end_fdate <- strftime(as.Date(end_date, origin="1970-01-01", tz="GMT"), format = "%Y%m%d", tz = "GMT", usetz = FALSE)
  f_list <- f_list_full[which(f_list_full_date >= start_fdate & f_list_full_date <= end_fdate)]
  # acq_date_list <- f_list_full_date[which(f_list_full_date >= start_fdate & f_list_full_date <= end_fdate)]
  
  # define reference EPSG
  r <- terra::rast(f_list[1])
  ref_epsg <- as.integer(terra::crs(r, describe=TRUE)$code)
  
  # import polygons
  p_in <- sf::st_read(input_polygon, stringsAsFactors = FALSE, quiet = TRUE)
  p_in$pid <- as.integer(1:nrow(p_in))
  
  # reproject polygons to target raster incidence angle
  if(ref_epsg != sf::st_crs(p_in)$epsg){
    p_in <- sf::st_transform(p_in, ref_epsg)
  }
  
  # clean workspace
  suppressWarnings(rm(list=c("f_list_full","f_list_full_date","start_fdate","end_fdate","ref_epsg","start_date","end_date","input_path","input_polygon")))
  invisible(gc())
  
  # define reference raster
  if(terra::res(r)[1] != sp_res){
    r_void <- terra::rast(extent=terra::ext(r), resolution=sp_res, crs=terra::crs(r))
  } else {
    r_void <- terra::rast(extent=terra::ext(r), resolution=terra::res(r)[1], crs=terra::crs(r))
  }
  # rasterize polygons
  p_boundary <- terra::vect(sf::st_cast(p_in, "MULTILINESTRING"))
  pl <- terra::rasterize(p_boundary, r_void, field=0, background=1, touches=TRUE)
  # pc <- terra::rasterize(p_in, r, field="pid", fun="count")
  r_cover <- terra::rasterize(p_in, r_void, cover=TRUE)
  r_cover <- ceiling(r_cover * pl)
  suppressWarnings(rm(list=c("pl","p_boundary")))
  invisible(gc())
  # p <- terra::rasterize(p_in, r_void, field="pid")
  p <- terra::rasterize(p_in, r_void, field="pid")
  p <- terra::ifel(r_cover == 1, p, NA)
  suppressWarnings(rm(list=c("r_cover","p_in")))
  invisible(gc())
  if(terra::res(p)[1] != terra::res(r)[1]){
    p <- terra::resample(p, r, method="near")
  }
  suppressWarnings(rm(list=c("r","r_void","sp_res")))
  invisible(gc())
  if(!(terra::inMemory(p))){
    p <- terra::toMemory(p)
  }
  
  # setup cycle
  ll <- length(f_list)
  out_df <- data.frame("pid"=integer(), "date"=character(), "epoch"=integer(), "B2"=double(), "B3"=double(), "B4"=double(), "B5"=double(), "B6"=double(), "B7"=double(), "B8"=double(), "B8A"=double(), "B11"=double(), "B12"=double())
  write.table(out_df, file=output_file, sep=";", dec=".", row.names=FALSE, col.names=TRUE, na="", quote=FALSE)
  invisible(gc())
  
  # initialize progress bar
  if(progress){
    pb <- utils::txtProgressBar(style=3)
  }
  
  # cycle on input files
  for(f in 1:ll){
    r <- terra::rast(f_list[f])
    
    # query spectral bands (raster-based approach)
    qut <- terra::zonal(r, p, fun="median", na.rm=TRUE)
    # rename ID field
    names(qut)[which(names(qut) %in% "layer")] <- "pid"
    
    # # query spectral bands (vector-based approach)
    # qut <- terra::extract(r, p_in, fun="median", touches=FALSE, small=FALSE, na.rm=TRUE)
    # # rename ID field
    # names(qut)[which(names(qut) %in% "ID")] <- "pid"
    
    # set acquisition date
    fname <- basename(f_list[f])
    date_dt_format <- paste(substr(fname, 8, 11), "-", substr(fname, 12, 13), "-", substr(fname, 14, 15), " ", substr(fname, 17, 18), ":", substr(fname, 19, 20), ":", substr(fname, 21, 22), sep="")
    epoch_value <- as.numeric(as.POSIXct(strftime(date_dt_format, format="%Y-%m-%d %H:%M:%S", tz="GMT", usetz = FALSE), tz = "UTC"))
    qut$date <- as.character(date_dt_format)
    qut$epoch <- as.integer(epoch_value)
    
    # # bind result to output data frame
    # out_df <- rbind(out_df, qut[,c(1,12,13,2,3,4,5,6,7,8,9,10,11)])
    
    # append results to output file
    write.table(qut[which(complete.cases(qut)),c(1,12,13,2,3,4,5,6,7,8,9,10,11)], file=output_file, sep=";", dec=".", row.names=FALSE, col.names=FALSE, na="", quote=FALSE, append=TRUE)
    
    suppressWarnings(rm(list=c("r","qut","fname","date_dt_format","epoch_value")))
    invisible(gc())
    
    if(progress){
      utils::setTxtProgressBar(pb, f/ll)
    }
  }
  
  if(progress){
    close(pb)
  }
  
  # # export result
  # write.table(out_df, file=output_file, sep=";", dec=".", row.names=FALSE, col.names=TRUE, na="", quote=FALSE)
  
  message(paste(c("Query ended at: "), Sys.time(), sep=""))
  message(paste("Elapsed time ", as.character(paste(as.integer(as.numeric(proc.time() - ptm)[3]/3600), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/3600 - as.integer(as.numeric(proc.time() - ptm)[3]/3600)) * 60)), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/60 - as.integer(as.numeric(proc.time() - ptm)[3]/60)) * 60)), sep="")), " hours", sep=""))
  
  # clean workspace and exit
  suppressWarnings(rm(list=ls()))
  invisible(gc())
  q("no")
  
}

#s2_select_measures_and_mean <- function()



