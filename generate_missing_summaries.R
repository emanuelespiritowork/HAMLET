#generate_summaries.R
source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/matchup.R")

#choose output summary folder
output_summary_folder = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Results/summary/"

#choose PRISMA level for matchups (L1 with ISOFIT or L2 standard)
#L1 with ISOFIT
# repo_root_path_prs <- "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/processing_L1" 
#or L2 standard
repo_root_path_prs = "//10.0.1.243/nr_data/3_rs_data/PRISMA/JDS/put_PRISMA_he5_and_S2_tif_here_L2"

generate_all_missing_summaries(output_summary_folder = output_summary_folder,
                               repo_root_path_prs = repo_root_path_prs)
