#main.R algorithm
#Author: Emanuele Spirito
#place: CNR-IREA

#load libraries 
source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/matchup.R")

#Parameter of selection:  
date_start = "2025-05-01"
date_end = "2025-05-10"
output_file = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results"
hypstar_loc = sf::st_read()
pivoted = T #data to be pivoted (true) or not (false)

#get hypstar data
hypstar_data <- get_hypstar_data(date_start = date_start, 
                                 date_end = date_end, 
                                 pivoted = pivoted)

#prepare hypstar for matchup
hypstar_for_matchup <- hypstar_select_measures_and_mean(dataset_res = hypstar_data,
                                                        input_pivoted = pivoted, 
                                                        output_pivoted = pivoted)

#s2_data 

s2_data <- get_s2_data(input_polygon = ?,
                       output_file = output_file,
                       start_date = date_start,
                       end_date = date_end) 
  
#s2_for_matchup


