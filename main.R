#main.R algorithm
#Author: Emanuele Spirito
#place: CNR-IREA

#load libraries 
source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/matchup.R")

#Parameter of selection:  
library(tidytable)
date_start <- "2025-05-01"
date_end <- "2025-05-10"
output_s2_file <- "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/s2_data.csv"
output_hypstar_file <- "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/hypstar_data.csv"
pivoted <- F #data to be pivoted (true) or not (false)
unpivoted_with_bands <- T #data to be unpivoted with bands (true) or not (false)
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
                                 date_end = date_end, 
                                 output_pivoted = pivoted,
                                 unpivoted_with_bands = unpivoted_with_bands)

#prepare hypstar for matchup
hypstar_angle_filtered <- hypstar_angle_filter_and_mean(dataset_res = hypstar_data,
                                                     input_pivoted = pivoted,
                                                     output_pivoted = pivoted,
                                                     unpivoted_with_bands = unpivoted_with_bands) 

write.table(hypstar_data, output_hypstar_file, quote = F, sep = ";",
            row.names = F)

#get s2 data 

s2_data <- get_s2_data(date_start = date_start,
                       date_end = date_end,
                       output_pivoted = pivoted,
                       unpivoted_with_bands = unpivoted_with_bands) 

#prepare s2 for matchup (optional)
s2_angle_filtered <- s2_angle_filter(dataset_res = s2_data,
                                     input_pivoted = pivoted,
                                     output_pivoted = pivoted,
                                     unpivoted_with_bands = unpivoted_with_bands)

write.table(s2_data, output_s2_file, quote = F, sep = ";",
            row.names = F)

matchup_result <- matchup_function(s2_data = s2_data,
                                   hypstar_data = hypstar_data,
                                   input_pivoted = pivoted,
                                   output_pivoted = pivoted,
                                   unpivoted_with_bands = unpivoted_with_bands)


############################
trasporre <- read.csv("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/s2_bands.csv", sep = ";")
aggiunta <- colnames(trasporre)[-1]
datatable <- data.table::data.table() %>%
  mutate(prima_colonna = aggiunta) %>%
  mutate(seconda_colonna = as.character(trasporre[-1])) %>%
  mutate(prima_colonna = as.numeric(gsub("X","",prima_colonna))) %>%
  rename(wvl = prima_colonna) %>%
  rename(bands = seconda_colonna)

write.table(datatable, "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/s2_bands_wavelenght/s2_bands.csv",
            sep = ",", quote = F, row.names = F)

