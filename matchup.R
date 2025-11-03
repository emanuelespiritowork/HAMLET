#Create match-up

#use res variable to make matchup with satellite image

select_measures <- function(dataset_res, #dataset pivoted of (date,time,VZA,VAA,PAA,B*)
                            minVZA, #minimum View Zenith Angle
                            maxVZA, #maximum View Zenith Angle
                            minVAA, #minimum View Azimuth Angle
                            maxVAA, #maximum View Azimuth Angle
                            delVZA, #vector deleting View Zenith Angle inside dataset
                            delVAA #vector deleting View Azimuth Angle inside datasets
                            ){
  library(tidytable)
  filtered_res <- dataset_res %>% filter(VZA > minVZA & VZA < maxVZA) %>%
    filter(VAA > minVAA & VAA < maxVAA) %>%
    filter(!(VZA %in% delVZA)) %>%
    filter(!(VAA %in% delVAA))
  
  matchup_selection_res <- filtered_res %>% 
    group_by(date,time,bands) %>%
    mutate(bands_mean = mean(bands)) %>%
    ungroup() %>%
    select(-bands) %>%
    rename(bands = bands_mean)
}
