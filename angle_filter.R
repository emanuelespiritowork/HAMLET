#angle_filter
#Author: Emanuele Spirito @ CNR-IREA
#Contributions: Gabriele Candiani, Lorenzo Parigi, Monica Pepe, Mirco Boschetti, Vittorio Brando
#Year: 2025
#Purpose: filter values in hypstar station and satellite data

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
