#main.R algorithm
#Author: Emanuele Spirito
#place: CNR-IREA

#load libraries 
source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/matchup.R")
#source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/F01_Spectral_Resampling_V02.R")

#Parameter of selection:  
# library(tidytable)
# date_start <- "2024-11-12"
# date_end <- "2024-11-14"
# date_start <- "2025-04-20"
# date_end <- "2025-04-30"
date_start <- "2025-04-23"
date_end <- "2025-04-25"
repo_root_path_prs <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/"
PRISMA_level_product <- "L2D"
processing_level_prs <- "coreg" #or 0_read
mirai::daemons(0)
mirai::daemons(4)

#minVZA = ?, #minimum View Zenith Angle
#maxVZA = ?, #maximum View Zenith Angle
#minVAA = ?, #minimum View Azimuth Angle
#maxVAA = ?, #maximum View Azimuth Angle
#minSZA = ?, #minimum Solar Zenith Angle
#maxSZA = ?, #maximum Solar Zenith Angle
#minSAA = ?, #minimum Solar Azimuth Angle
#maxSAA = ?, #maximum Solar Azimuth Angle
#delVZA = c(""), #vector deleting View Zenith Angle inside dataset
#delSZA = c(""), #vector deleting Solar Zenith Angle inside dataset
#delVAA = c(""), #vector deleting View Azimuth Angle inside datasets


#get hypstar data
hypstar_data <- get_hypstar_data(date_start = date_start, 
                                 date_end = date_end)

#get s2 data 
s2_data <- get_s2_data(date_start = date_start,
                       date_end = date_end)

#matchup to S2

matchup_result_to_S2 <- matchup_function(sat_data = s2_data,
                 hypstar_data = hypstar_data,
                 time_or_sun_angle = T,
                 just_closest_time_or_sun_angle = T,
                 just_closest_view_angle = F,
                 range_view_angle = 360,
                 delVZA = c(0),
                 weighted_mean = F,
                 print_polar = T)

plot <- plot_ribbon(matchup_result_to_S2)
plot <- plot_rainbow(matchup_result_to_S2)
plot <- plot_relative_errors(matchup_result_to_S2)

#get PRISMA data 

PRISMA_data <- get_PRISMA_data(date_start = date_start,
                               date_end = date_end,
                               repo_root_path = repo_root_path_prs,
                               processing_level = processing_level_prs,
                               PRISMA_level_product = PRISMA_level_product) 

#matchup to PRISMA
matchup_result_to_PRISMA <- matchup_function(sat_data = PRISMA_data,
                                             hypstar_data = hypstar_data,
                                             time_or_sun_angle = T,
                                             just_closest_time_or_sun_angle = T,
                                             just_closest_view_angle = F,
                                             range_view_angle = 30,
                                             delVZA = c(0),
                                             weighted_mean = F,
                                             print_polar = T)

plot <- plot_ribbon(matchup_result_to_PRISMA)
plot <- plot_rainbow(matchup_result_to_PRISMA)
plot <- plot_relative_errors(matchup_result_to_PRISMA)

matchup_result_to_PRISMA <- matchup_function(sat_data = PRISMA_data,
                                             hypstar_data = hypstar_data,
                                             time_or_sun_angle = T,
                                             just_closest_time_or_sun_angle = F,
                                             range_time_or_sun_angle = 2,
                                             just_closest_view_angle = F,
                                             range_view_angle = 360,
                                             delVZA = c(0),
                                             weighted_mean = F,
                                             print_polar = T)

plot <- plot_ribbon(matchup_result_to_PRISMA)
plot <- plot_rainbow(matchup_result_to_PRISMA)
plot <- plot_relative_errors(matchup_result_to_PRISMA)

matchup_result_to_PRISMA <- matchup_function(sat_data = PRISMA_data,
                                             hypstar_data = hypstar_data,
                                                     time_or_sun_angle = T,
                                                     just_closest_time_or_sun_angle = T,
                                                     range_time_or_sun_angle = NULL,
                                                     just_closest_view_angle = F,
                                                     range_view_angle = 30,
                                                     delVZA = c(0),
                                             weighted_mean = F,
                                             print_polar = T)

plot <- plot_ribbon(matchup_result_to_PRISMA)
plot <- plot_rainbow(matchup_result_to_PRISMA)
plot <- plot_relative_errors(matchup_result_to_PRISMA)

#get EnMAP data

EnMAP_data <- get_EnMAP_data(date_start = date_start,
                             date_end = date_end)

#matchup to EnMAP
matchup_result_to_EnMAP <- matchup_function(sat_data = EnMAP_data,
                                            hypstar_data = hypstar_data,
                                                     time_or_sun_angle = T,
                                                     just_closest_time_or_sun_angle = F,
                                                     range_time_or_sun_angle = 1,
                                                     just_closest_view_angle = F,
                                                     range_view_angle = 360,
                                                     delVZA = c(0),
                                            print_polar = T)

plot <- plot_ribbon(matchup_result_to_EnMAP)
plot <- plot_rainbow(matchup_result_to_EnMAP)
plot <- plot_relative_errors(matchup_result_to_EnMAP)

matchup_result_to_EnMAP <- matchup_function(sat_data = EnMAP_data,
                                            hypstar_data = hypstar_data,
                                                     time_or_sun_angle = T,
                                                     just_closest_time_or_sun_angle = T,
                                                     range_time_or_sun_angle = NULL,
                                                     just_closest_view_angle = F,
                                                     range_view_angle = 30,
                                                     delVZA = c(0))

plot <- plot_ribbon(matchup_result_to_EnMAP)
plot <- plot_rainbow(matchup_result_to_EnMAP)
plot <- plot_relative_errors(matchup_result_to_EnMAP)

# ############################
# write.table(matchup_result_to_PRISMA,"//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/matchupPRISMA_1hour_range.csv",
#             sep = ";", quote = F, row.names = F, append = F)
# 
# write.table(matchup_result_to_EnMAP,"//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/csv/matchupEnMAP_1hour_range.csv",
#             sep = ";", quote = F, row.names = F, append = F)
# 
# ############################
# trasporre <- read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/s2_bands.csv", sep = ";")
# aggiunta <- colnames(trasporre)[-1]
# datatable <- data.table::data.table() %>%
#   mutate(prima_colonna = aggiunta) %>%
#   mutate(seconda_colonna = as.character(trasporre[-1])) %>%
#   mutate(prima_colonna = as.numeric(gsub("X","",prima_colonna))) %>%
#   rename(wvl = prima_colonna) %>%
#   rename(bands = seconda_colonna)
# 
# write.table(datatable, "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/s2_bands.csv",
#             sep = ",", quote = F, row.names = F)
# ############################
# s2a <- read.table("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/S2C_SRF_10b.csv", sep = ",")
# names(s2a) <- c("hyp_wvl","B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12")
# s2a <- s2a[-1,]
# write.table(s2a,"//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/S2C_SRF_10b_bands.csv", sep = ",",
#                    quote = F, row.names = F)
# 
# ############################
# s2c <- read.table("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/config_data/SRF/S2C_SRF_13b.csv",
#                   sep = ",") %>% janitor::row_to_names(row_number = 1)
# names(s2c) <- c("band","B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B10","B11","B12")
# s2c_fill <- data.frame(1:299,0,0,0,0,0,0,0,0,0,0,0,0,0)
# names(s2c_fill) <- c("band","B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B10","B11","B12")
# s2c <- s2c_fill %>% bind_rows(s2c) 
# s2c <- s2c %>% mutate(band = as.numeric(band)) 
# write.table(s2c, "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/config_data/SRF/S2C_SRF_13b.csv",
#             sep = ",", quote = F, row.names = F)
# 
# ##############################
# Legenda
# Enmap color = "red"
# Hypstar color = "blue"
# PRISMA color = "green"

ncin$var

function(var){
  return(c(var$name,var$longname))
}
