library(dplyr)
library(glue)

source("./help_simul.R")

D = 1000  # Number of simulations
m = 400   # Number of clusters per each simulation

dir.create(glue("data/m{m}"), recursive = TRUE)

set.seed(1)


start_time = Sys.time()
progress_marks <- floor(seq(0.1, 1, by = 0.1) * D)
mark_idx <- 1

cat("\n")
for(i in 1:D){
  
  ## Show progress bar for every 10%
  if (mark_idx <= length(progress_marks) && i == progress_marks[mark_idx]) {
    pct <- mark_idx * 10
    et <- round(difftime(Sys.time(), start_time, units = "secs"), 1)
    bar <- paste0("[", strrep("*", mark_idx), strrep("-", 10 - mark_idx), "]")
    cat(glue("   {bar} {pct}% | ET: {et}s"), "\n")
    mark_idx <- mark_idx + 1
  }
  
  data = data.sim(m)
  saveRDS(data, glue("data/m{m}/data_id{i}.rds"))
  
}


