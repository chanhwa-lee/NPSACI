time = proc.time()

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(tidyverse)
library(optparse)

###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-M", "--M"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters for estimand simulation")),
  make_option(c("-r", "--r"), action = "store", default = NA, type = "integer",
              help = paste0("Number of binary vector sampling for outcome reg computation")),
  make_option(c("-n", "--n"), action = "store", default = NA, type = "integer",
              help = paste0("N distribution scenario.
                            1: N ~ Unif{3,...,5}
                            2: N ~ Unif{5,...,20}
                            3: N ~ Unif{20,...,50}
                            4: N ~ Unif{50,...,100}
                            5: N ~ CholeraVaccineData"))
)

opt = parse_args(OptionParser(option_list = option_list))

M = opt$M              # Number of clusters for target estimand computation
r = opt$r              # Number of binary vector sampling
ndist = opt$n          # N distribution scenario

taus = 0.2
thetas = c(0.3,0.4,0.5,0.6,0.7)
theta0 = 0.5

###----------- Source help functions  ---------------###
source("help_simul.R")

###----------- Print simulation setting  ---------------###
cat("\n[Simulation setting]","\n\n")
cat(paste0("- M: ", M),"\n\n")
cat(paste0("- taus: ", paste(signif(taus, 4), collapse = ", ")),"\n\n")
cat(paste0("- thetas: ", paste(signif(thetas, 4), collapse = ", ")),"\n\n")
cat(paste0("- r: ", r),"\n\n")

###----------- Compute estimands  ---------------###
estimands = estimands.sim(M, "TypeB", taus, thetas, theta0, r, ndist)
print(estimands)

###-------- Save simulated estimand list as Rdata --------###

## Save output
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
saveRDS(estimands, file = glue("estimand/Ndist{ndist}/estimand_id{task_ID}.rds"))

## Stop timer and report total run time
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))