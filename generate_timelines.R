#generate_timeline.R
source("//10.0.1.243/nr_working/emanuele/Progetto_JSIT/JSIT-hypstar/matchup.R")

#insert timeline repo of dates
timeline_repo = "//10.0.1.243/nr_working/emanuele/Progetto_JSIT/Timelines/csv"

#generate timeline
timeline <- create_timeline(timeline_repo = timeline_repo)

timeline
