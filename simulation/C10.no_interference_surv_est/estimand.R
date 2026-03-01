time = proc.time()

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.5")
library(dplyr)
library(tidyverse)
library(doMC)
library(foreach)
library(optparse)

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

taus = 0.2
thetas = c(0, 0.5, 1.0)
theta0 = 0

## Number of cpus (cores) for parallel computing
n.cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
n.cpus <- ifelse(is.na(n.cpus), 3, n.cpus)
registerDoMC(cores = n.cpus) 

###----------- Source help functions  ---------------###
source("help_simul.R")

###----------- Print simulation setting  ---------------###
cat("\n[Simulation setting]","\n\n")
cat(paste0("- M: ", M),"\n\n")
cat(paste0("- taus: ", paste(signif(taus, 4), collapse = ", ")),"\n\n")
cat(paste0("- thetas: ", paste(signif(thetas, 4), collapse = ", ")),"\n\n")
cat(paste0("- r: ", r),"\n\n")

###----------- Compute estimands  ---------------###
estimands = estimands.sim(M, "TypeB", taus, thetas, theta0, r)
print(estimands)

###-------- Save simulated estimand list as Rdata --------###

## Save output
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
saveRDS(estimands, file = glue("estimand/TypeB_alpha_0_1/estimand_id{task_ID}.rds"))

## Stop timer and report total run time
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))