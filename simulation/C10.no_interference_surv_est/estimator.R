time = proc.time()

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.5")
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
              help = paste0("Number of clusters per each simulation")),
  make_option(c("-r", "--r"), action = "store", default = NA, type = "integer",
              help = paste0("Number of binary vector sampling for outcome reg computation"))
)

opt = parse_args(OptionParser(option_list = option_list))

m = opt$m              # Number of clusters per simulation
r = opt$r              # Number of binary vector sampling

## Wrap TASK_ID with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

taus = 0.2
thetas = c(0, 0.5, 1)
theta0 = 0
  
## Number of cpus (cores) for parallel computing
n.cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
n.cpus <- ifelse(is.na(n.cpus), 3, n.cpus)
registerDoMC(cores = n.cpus) 


###----------- Source help functions  ---------------###
code.dir = "~/research/NPSACI/code"

## Help functions for simulated data
source("../../help_simul.R")

## Help functions for estimator main functions
source(paste0(code.dir,"/help_util.R"))

## Help functions for nuisance functions estimation method
source(paste0(code.dir,"/help_nuis_est.R"))

###----------- Import dataset  ---------------###
data = data.sim(m)

###-------- Fit estimator --------###

## Covariate names
X.names = colnames(data %>% select(-c(id, Y, D, A)))

## Number of sample splints
K = 2

result = estimator(data = data,
                   X.T.names = X.names, X.C.names = X.names, X.A.names = X.names,
                   policy = "TypeB",
                   taus = taus, thetas = thetas, theta0 = theta0, 
                   K = K, r = r, bounding = TRUE,
                   Tfit.method = "nonpara", Cfit.method = "nonpara", Afit.method = "nonpara")

###-------- Save simulated estimator list as Rdata --------###

## Save output
saveRDS(result$result, file = paste0("Rdata/estimate_id", task_ID,".rds"))

## Stop timer and report total run time
script.time = proc.time() - time
cat(glue("\n--- Total run time: {round(script.time[3], 1)} seconds --- \n"))
