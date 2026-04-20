# HEADER --------------------------------------------
#
# Author: 		Lorenzo Parigi
# Email:  		parigi.l@irea.cnr.it
#
# Date: 		15/03/2024
#
# Script Name:	F01_Spectral_Resampling_V02.R
#
# Version:      2
# 
#
# Script Description: The script allows to resample a dataframe or a raster based on a Spectral Response function or using the Full Width Half Maximum of a sensor.
#
#
####Instructions####
# 
# The first step is to call the function SRF_adapter(origin_band_center, mode, sensor). The origin_band_center must be a numeric vector with the central wavelength of the original data, the mode can be SRF or FWHM based on the tipe of processing required, the sensor is the code of the target sensor. The output of the function is a matrix and it must be saved.
# Example: SFR_S2A <- SRF_adapter(c(400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500), mode = "SRF", sensor = "S2A")
# The second step is to call raster_resampler(origin_raster, SRF) or dataframe_resampler(total_df, band_variables, SRF). In the first case the raster and the SRF_matrix produced by the SFR_adapter function are enough. The dataframe_resampler() needs the dataset, and the names of the columns corresponding to the bands that have to be resampled.
# 
# It is possible to input a custom SRF or FWHM configuration, in this case you need to write "custom" as sensor and input a dataframe in the SRF_DB or FWHM_DB function parameter.
# Example: SFR_S2A<- SRF_adapter(c(400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500), mode = "SFR", sensor = "custom", SRF_DB = custom_S2A_DF)
# 
#
#



####Functions####
#####Resampling function#####
resampling_bands <- function(central_band, band_list, value_list) {
    weights_list <- 1/(abs(central_band-band_list+0.00000001))^2
    
    value_end <- sum(value_list*weights_list)/sum(weights_list)
    
    return(value_end)
}

#####Wrapper#####
resampling_wrapper_function <- function(goal_bands, band_list, value_list) {
    result <- lapply(goal_bands, function(central_band) {
        resampling_bands(central_band, band_list, value_list)
    }) |> unlist()
    
    return(result)
}

#####SFR producer#####
SRF_calculator <- function(SRF_DB, band_centers) {
    SRF <- SRF_DB |>
        tidytable::reframe(tidytable::across(.cols = -c(band),
                                             .fns = ~ resampling_wrapper_function(band_centers, band, .x))) |>
        tidytable::mutate(tidytable::across(tidytable::everything(),
                                            .fns = ~ .x/sum(.x))) |>
        as.matrix()
    
    colnames(SRF) <- setdiff(names(SRF_DB), "band")
    
    return(SRF)
}

#####FWHM calculator#####
FWHM_calculator <- function(FWHM_DB, band_centers) {
    
    result <- lapply(1:length(FWHM_DB[["center"]]), function(x) {
        
        dnorm(band_centers, FWHM_DB[["center"]][x], FWHM_DB[["fwhm"]][x]/(2*sqrt(2*log(2))))
    }) |>
        data.table::as.data.table() |>
        tidytable::mutate(tidytable::across(tidytable::everything(),
                                            .fns = ~ .x/sum(.x))) |>
        as.matrix()
    
    colnames(result) <- FWHM_DB[["Name"]]
    
    return(result)
}

#####SRF Matrix Generator#####
SRF_adapter <- function(origin_band_centers = NULL,
                        mode = "SRF", #SRF or FWHM
                        sensor = "S2A", #Sensor name or custom
                        ref_folder = here("config_data"),
                        band_as_title = FALSE, 
                        SRF_DB = NULL, #if the sensor is specified, these param are not needed
                        FWHM_DB = NULL) {
    if (mode == "SRF") {
        if (sensor != "custom") {
            SRF_DB <- list.files(paste0(ref_folder, "/SRF"), full.names = TRUE) 
            SRF_DB <- SRF_DB[grepl(sensor, SRF_DB)] |>
                tidytable::fread()
        }
        
        SRF <- SRF_calculator(SRF_DB, origin_band_centers)
        
    } else if (mode == "FWHM") {
        if (sensor != "custom") {
            FWHM_DB <- list.files(paste0(ref_folder, "/FWHM"), full.names = TRUE) 
            FWHM_DB <- FWHM_DB[grepl(sensor, FWHM_DB)] |>
                tidytable::fread()
        }
        
        SRF <- FWHM_calculator(FWHM_DB, origin_band_centers)
        
    } else {
        stop("mode must be SRF or FHWM")
    }
    
    
    return(SRF)
    
}

#####Raster wrapper#####
raster_resampler <- function(origin_raster, SRF) {
    
    raster_dimensions <- dim(origin_raster) 
    raster_dimensions <-raster_dimensions[c(2,1)] 
    raster_dimensions <- c(raster_dimensions, ncol(SRF))
    
    raster_extent <- terra::ext(origin_raster)
    raster_crs <- terra::crs(origin_raster)
    
    final_raster <- origin_raster |>
        as.matrix()%*%SRF 
    
    final_raster <- final_raster |>
      array(raster_dimensions) |>
      terra::rast() |> 
      t() |>
      terra::`crs<-`(raster_crs) |>
      terra::`ext<-`(raster_extent)
    
    names(final_raster) <- colnames(SRF)
    
    return(final_raster)
    
}


#####Resamler#####
dataframe_resampler <- function(total_df, band_variables, SRF) {
    
    variable_df <- total_df |>
        tidytable::select(tidytable::any_of(setdiff(names(total_df), band_variables)))
    bands_df <- total_df |> 
        tidytable::select(tidytable::any_of(band_variables)) |>
        as.matrix()
    
    result <- bands_df%*%SRF |>
        tidytable::data.table() |>
        `names<-`(colnames(SRF)) 
    result <- cbind(variable_df, result)
    
    return(result)
    
}

####Versions####
#Improved the raster resampling


####TO-DO####
# Add a dataframe with all the sensor and the path to the configuration files
# Point to the dataframe when looking for a sensor
# Add a function that output the DF with the name of the sensors
# Comment the script