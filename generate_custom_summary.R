#generate custom summary
source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/matchup.R")
mirai::daemons(0)
mirai::daemons(8)
#choose output summary folder
output_summary_folder = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/summary/"

#choose PRISMA level for matchups (L1 with ISOFIT or L2 standard)
#L1 with ISOFIT
# repo_root_path_prs <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/processing_L1" 
#or L2 standard
repo_root_path_prs = "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/put_PRISMA_he5_and_S2_tif_here_L2"

#summarytype can be bestconfig or average
summary_type = "bestconfig"

#hypstar level can be L2A or L2B
hypstar_level = "L2B"

create_summary(date_start = "2025-03-01",
               date_end = "2025-06-30",
               repo_root_path_prs = repo_root_path_prs,
               output_summary_folder = output_summary_folder,
               hypstar_level = hypstar_level,
               summary_type = summary_type)
