time = proc.time()

###------------------- Load libraries ----------------------###
library(dplyr)
library(glue)
library(survival)
library(doMC)
library(randomForestSRC)
library(reshape2)
library(dbarts)
library(optparse)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(reshape2::melt)


###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-m", "--m"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters per each simulation"))
)

opt = parse_args(OptionParser(option_list = option_list))

m = opt$m              # Number of clusters per simulation

## Wrap TASK_ID with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## TypeB policy parameters
taus = 0.2

###----------- Source help functions  ---------------###
code.dir = "~/research/NPSACI/code"

## Help functions for simulated data
source("../../help_simul.R")
source("../../RSF_learners.R")

###----------- Import dataset  ---------------###
data = data.sim(m) %>%
  dplyr::rename(Time = Y, Event = D, Treatment = A) %>%
  dplyr::select(-id)

###-------- Fit estimator --------###
df_results <- get_all_ATE_results(data, time.interest=taus)

print("--------------------------------------------------")
print(" FINAL ESTIMATED ATE RESULTS (95% CI) ")
print("--------------------------------------------------")
print(df_results)


###-------- Save simulated estimator list as Rdata --------###

## Save output
saveRDS(df_results, file = glue("Rdata/estimate_id{task_ID}.rds"))

## Stop timer and report total run time
script.time = proc.time() - time
cat(glue("\n--- Total run time: {round(script.time[3], 1)} seconds --- \n\n"))
