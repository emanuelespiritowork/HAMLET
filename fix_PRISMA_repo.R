#copy all_angles_file_to_PRISMA_L2

find_all_angles_file_L1 <- list.files(path = "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/processing_L1",
                                      pattern = "all_angles_file.csv",
                                      recursive = T,
                                      full.names = T)

from_file <- paste0(list.dirs(path = "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/processing_L1",
          recursive = F,
          full.names = T),"/PRISTAR_processing/0_read/ATCOR/all_angles_file.csv")

to_file <- gsub( "PRS_L2C_STD_2026", "PRS_L2D_STD_2026",
  gsub( "PRS_L2C_STD_2025", "PRS_L2D_STD_2025",
  gsub( "PRS_L1_STD_OFFL_", "PRS_L2C_STD_",
  gsub("processing_L1","put_PRISMA_he5_and_S2_tif_here_L2",from_file) )))

file.copy(from_file, to_file, overwrite = T)
