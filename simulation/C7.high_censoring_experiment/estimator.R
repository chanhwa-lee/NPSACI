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
              help = paste0("Number of clusters per each simulation")),
  make_option(c("-r", "--r"), action = "store", default = NA, type = "integer",
              help = paste0("Number of binary vector sampling for outcome reg computation")),
  make_option(c("-c", "--c"), action = "store", default = NA, type = "double",
              help = paste0("Target censoring rate in simulated data"))
)

opt = parse_args(OptionParser(option_list = option_list))

m = opt$m              # Number of clusters per simulation
r = opt$r              # Number of binary vector sampling
censoring_rate = round(opt$c, 1) # Target censoring rate in simulated data

## Table of censoring_factor and censoring_rate ##
# censoring_factor | censoring_rate
# 0.91 | 0.9
# 0.81 | 0.8
# 0.67 | 0.7
# 0.51 | 0.6
# 0.28 | 0.5
# 0.00 | 0.4

# Rounding to 1 decimal place to ensure matching
censoring_factor = case_when(
  censoring_rate == 0.9 ~ 0.91,
  censoring_rate == 0.8 ~ 0.81,
  censoring_rate == 0.7 ~ 0.67,
  censoring_rate == 0.6 ~ 0.51,
  censoring_rate == 0.5 ~ 0.28,
  censoring_rate == 0.4 ~ 0.00,
  TRUE        ~ as.numeric(NA)
)

cat(glue("\n--- Target Censoring rate: {censoring_rate} || Censoring factor: {censoring_factor} --- \n\n"))

if (is.na(censoring_factor)) {
  stop(paste0("Censoring rate '", censoring_rate, "' is not specified correctly."))
}

## Wrap TASK_ID with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## TypeB policy parameters
taus = c(0.2, 0.4)
thetas = 0.3 + 0.01*0:30
theta0 = 0.3 + 0.01*15 ## Prevent potential floating issue 0.3+0.05*15 != 0.45

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
data = data.sim(m, censoring_factor = censoring_factor)
cat(glue("\n--- Censoring rate: {round(1-mean(data$D),3)*100}% --- \n\n"))
cat(glue("\n--- Desired censoring rate: {censoring_rate*100}% --- \n\n"))

###-------- Fit estimator --------###

## Covariate names
X.names = colnames(data %>% dplyr::select(-c(id, Y, D, A)))

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
saveRDS(result$result, file = glue("Rdata/estimate_id{task_ID}.rds"))

## Stop timer and report total run time
script.time = proc.time() - time
cat(glue("\n--- Total run time: {round(script.time[3], 1)} seconds --- \n\n"))
