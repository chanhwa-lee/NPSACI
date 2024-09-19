time = proc.time()

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(tidyverse)
library(survival)

suppressPackageStartupMessages(require(optparse))

## Help functions for simulated data
source("help_simul.R")


###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-M", "--M"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters for estimand simulation")),
  make_option(c("-r", "--r"), action = "store", default = NA, type = "integer",
              help = paste0("Number of binary vector sampling for outcome reg computation"))
)

opt = parse_args(OptionParser(option_list = option_list))

M = opt$M              # Number of clusters for target estimand computation
r = opt$r              # Number of binary vector sampling
taus = 0.1*(1:5)
thetas = 0.1*(0:10)
theta0 = 0

print("[Simulation setting]")
print(paste0("M: ", M))
print(paste0("taus: ", paste(signif(taus, 4), collapse = ", ")))
print(paste0("thetas: ", paste(signif(thetas, 4), collapse = ", ")))
print(paste0("r: ", r))

estimands = estimands.sim(M, taus, thetas, theta0, r)

print(estimands)

###-------- Save simulated estimand list as Rdata --------###
## Wrap with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Save output
saveRDS(estimands, file = paste0("estimand/estimand_id", task_ID,".rds"))


###------------------------------------------------------------###

## Stop timer and report total run time

script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))