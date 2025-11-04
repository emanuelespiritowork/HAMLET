rm(list=ls())

#### 0)Library ####

library(here)
library(ncdf4)
library(ggplot2)
library(dplyr)
library(FieldSpectroscopyDP)#remotes::install_github("tommasojulitta/FieldSpectroscopyDP")
library(FieldSpectroscopyCC)#remotes::install_github("tommasojulitta/FieldSpectroscopyCC")
library(tidytable)

#### 0)Function ####
smoothing_custom_function <- function(DB, 
                                      band_centers, #centri banda del sensore
                                      DF=60, #gradi di libertà
                                      ID_col = 1:5, #colonne del DB che non sono bande 
                                      start=c(580,1340), #inizio regione da escludere
                                      end=c(590,1430) #fine regione da escludere
){
  
  exclude_bands<-data.frame(start=start,end=end)
  
  res_excluderegions<-apply(DB %>% select(-ID_col),
                            1,function(x,y) {ExcludeSpectralRegions(x,y,SpectralRegion=exclude_bands)},x=band_centers)
  
  res_splineinterp<-apply(res_excluderegions,2,function(x,y) {SplineSmoothGapfilling(x,y,df=DF)},x=band_centers)
  res_splineinterp <- as.data.frame(res_splineinterp)
  res_splineinterp <- round(res_splineinterp,0)
  res_splineinterp_t <- as.data.frame(t(res_splineinterp))
  res_splineinterp_t <- cbind(DB %>% select(ID_col),res_splineinterp_t)
  colnames(res_splineinterp_t) <- colnames(DB)
  
  return(res_splineinterp_t)
}

#### 1) open files ####

# load the `ncdf4` and the `CFtime` packages 
home_path <- "//10.0.1.243/projects/2023_HyperNet-POP/3-DATA/processed"

# load HYPSTAR bands
bands_HYP <- read.csv("//10.0.1.243/projects/2023_HyperNet-POP/6-script/bands_HYP.csv")

##### 1.2) Parameter of selection:  start and end date ####
date_start = "2025-05-01"
date_end = "2025-05-10"

# creation of folder list
folder_list <- data.table::data.table(folder_name = list.dirs(home_path, recursive = FALSE, full.names = FALSE)) %>%
  mutate(date = substr(folder_name, 4, 11) %>% as.Date(., format = "%Y%m%d"),
         time = substr(folder_name, 12, 20)) %>%
  filter(date >= date_start,
         date <= date_end)

#### 2) create file of HYPSTAR reflectance ####

res <- lapply(1:nrow(folder_list), function(row_n){
    folder_name <- folder_list$folder_name[row_n]
    date <- folder_list$date[row_n] 
    time <- folder_list$time[row_n] 
    
    file_list <- list.files(paste0(home_path, "/", folder_name), pattern = "L2A_REF", 
                            full.names = TRUE)
    
    if (length(file_list) == 0) {
        return(NULL)
    }
    
    ncin <- nc_open(file_list)
    
    hyp_ref <-ncvar_get(ncin, "reflectance")
    hyp_wvl  <-ncvar_get(ncin, "wavelength")
    hyp_vaa <-ncvar_get(ncin, "viewing_azimuth_angle")
    hyp_paa   <-ncvar_get(ncin, "pointing_azimuth_angle")
    hyp_vza  <-ncvar_get(ncin, "viewing_zenith_angle")
    hyp_sza <- ncvar_get(ncin, "solar_zenith_angle")
    hyp_saa <- ncvar_get(ncin, "solar_azimuth_angle")
    
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

#### 3) interpolate spectra ####

res_smoothed <- smoothing_custom_function(DB = res,
                                          band_centers = bands_HYP$wvl, #centri banda del sensore
                                          DF=60,
                                          ID_col = 1:5, #colonne del DB che non sono bande 
                                          start=c(580,1340), #inizio regione da escludere
                                          end=c(590,1430))

#### 4) Manipulate data: create a long file spectra ####

# long for interpolated
res2int_long <- tidyr::pivot_longer(res_smoothed, cols = starts_with("B"),
                          names_to = "bands",
                          values_to = "Ref",
                          values_drop_na = TRUE) 
res2int_long <-left_join ( res2int_long,bands_HYP) %>%
  group_by(date, VZA, VAA, PAA, SZA, SAA, bands, wvl) %>%
  summarise(Ref = mean(Ref)) %>% ungroup()


# long for raw HYPSTAR data

res2_long <- tidyr::pivot_longer(res, cols = starts_with("B"),
                           names_to = "bands",
                           values_to = "Ref",
                           values_drop_na = TRUE) 
res2_long <-left_join ( res2_long,bands_HYP)

##### 4.1) manipulate long file ####

res2_long <-res2_long %>% mutate(month = substring(date, 6,7),
                                 year =  substring(date, 1,4),
                                 day = as.numeric(substring(date, 9,10)),
                                 week = ifelse(day < 14 , "1-14", "15-30"))

res2int_long <-res2int_long %>% mutate(month = strftime(date, format = "%m") %>% as.integer(),
                                       day = strftime(date, format = "%d") %>% as.integer(),
                                       week = strftime(date, format = "%V") %>% as.integer(),
                                       year = strftime(date, format = "%Y") %>% as.integer(),
                                       year = if_else(week == 1 & year == 2024, year+1, year),
                                       year_week = paste0(year, "_Week_", week)) %>%
  group_by(year_week) %>%
  mutate(min_date = min(date)) 

##### 4.2) check for anomalies ####

#anomalous_date <- res2int_long%>% filter(VZA == 30, VAA == 293, Ref >4500)

#anomalous_date <- unique(anomalous_date$date)

##### 4.3) Create file of weekly average ####

summarised_res_int <- res2int_long %>%
  group_by(VZA, VAA, PAA, bands, wvl, year_week, min_date) %>%
  summarise(Ref = mean(Ref)) %>% ungroup()

#### 5) Plotting ######

# setting: # define viewing geometry of interest

filter_VZA = 20
filter_VAA = 278

# 5.1) plot of spectra and weekly mean grouped by year_week #######

# order for facet

# "2024_Week_44","2024_Week_45","2024_Week_46","2024_Week_48","2024_Week_49",
# "2024_Week_50","2024_Week_51","2024_Week_52","2025_Week_1","2025_Week_2",
# "2025_Week_3","2025_Week_4","2025_Week_5",

neworder <- c("2025_Week_6","2025_Week_7","2025_Week_8","2025_Week_9","2025_Week_10", "2025_Week_11", "2025_Week_12", "2025_Week_13", "2025_Week_14",
              "2025_Week_15", "2025_Week_16", "2025_Week_17", "2025_Week_18", "2025_Week_19", "2025_Week_20", "2025_Week_21", "2025_Week_22", "2025_Week_23")

neworder <- c("2025_Week_18", "2025_Week_19" ,"2025_Week_20", "2025_Week_21" ,"2025_Week_22","2025_Week_23")

ggplot(data = summarised_res_int %>% filter(VZA==filter_VZA,
                                            VAA==filter_VAA), 
       aes(x = wvl, y = Ref)) +
  geom_point(data = res2int_long %>% filter(VZA==filter_VZA,
                                            VAA==filter_VAA), 
             color = "grey", alpha = 0.8)+
  geom_point(color = "red") +
  ylim(0,6000)+
  facet_wrap(~factor(year_week, levels =neworder ))+
  theme_light()+
  geom_label(data = res2int_long,
             aes(x = 400, y = 5000, label = min_date),
             hjust = 0, vjust = 1, inherit.aes = FALSE,
             fill = "white", alpha = 0.7, size = 3)

annotazioni <- df |>
  dplyr::group_by(year_week) |>
  dplyr::summarise(min_date = unique(min_date),
                   x = min(wvl, na.rm = TRUE),  # Puoi aggiustare per spostare il box
                   y = max(Ref, na.rm = TRUE))  # Idem, per metterlo in alto ## 5.2) plot of weekly mean on same graph #######

ggplot(data = summarised_res_int %>% 
         filter(VZA == filter_VZA,
                VAA == filter_VAA) %>%
         mutate(wvl =round(wvl, 3)) %>%
         arrange(wvl) %>%
         #select(-c(wvl)) %>%
         select(-c(bands)), 
       aes(x = wvl, y = Ref, color=as.factor(min_date))) +
  geom_point() + 
  # scale_fill_gradient(low="blue", high="red")+
  scale_color_gradientn(colours = rainbow(10))


## 5.3) plot single measurements on a range of date  #######

# define viewing geometry of interest
filter_VZA = 30
filter_VAA = 293


ggplot(data = res2int_long %>% filter(date > "2025-05-05", VZA==filter_VZA,
                                      #VAA==filter_VAA
                                      ) ,aes(x = wvl, y = Ref, color = as.factor(date))) +
  geom_point()+
  facet_wrap(~factor(VAA))+
  #facet_grid(factor(VAA)~as.factor(date))+
  theme_light()+
#  scale_color_gradientn(colours = rainbow(5))+
  scale_color_brewer(palette = "Spectral")

+
 # ylim(4000,6000)+
  #xlim(700,1200)

## 5.4) plot temporal trend of NIR @ 900 and NDVI  #######
# select data range
min_data <- "2025-05-01"
max_data <-  "2025-05-15" 
# select VAA
VZA_value <-   c(10,20,30)
# select VZA
VAA_value <- unique(res2int_long$VAA)
  Selected_band <- "B1128"
  


ggplot(data = res2int_long %>% filter(date > min_data,date < max_data,
                                      VZA %in% VZA_value,
                                      VZA %in% VZA_value,
                                       bands == Selected_band),
       aes(x = date, 
           y = Ref, 
           color = as.factor(date), 
           shape =as.factor(VZA))) +
  geom_point(aes(size=5))+
  geom_vline(xintercept = as.Date("2025-05-09"))+
  #ylim(4000,6000)+
  facet_grid(factor(VZA)~factor(VAA))
  #facet_wrap(~factor(VAA))+
  theme_light()
  
  

+
  scale_color_brewer(palette = "Spectral")



#### 6) Saving file ######

# define View and Zenith angle of interest
# VZA "0  5 10 20 30"
# VAA "83  98 113 263 278 293"

filter_VZA = 20
filter_VAA = 278

# long to wide

summarised_res_int_wide <- summarised_res_int %>% 
  filter(VZA == filter_VZA,
         VAA == filter_VAA) %>%
  mutate(wvl =round(wvl, 3)) %>%
  arrange(wvl) %>%
  #select(-c(wvl)) %>%
  select(-c(bands)) %>%
  pivot_wider(names_from = "wvl", values_from = "Ref")
write.csv(summarised_res_int_wide, "//10.0.1.243/projects/2023_HyperNet-POP/1_ACTIVITIES/TRAITS_RETRIEVAL/DB/VALIDATION/hypstar_week_AVG_241101_250309.csv")

