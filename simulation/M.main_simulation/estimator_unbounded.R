time = proc.time()

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
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
  make_option(c("-s", "--s"), action = "store", default = NA, type = "integer",
              help = paste0("Scenario indicating which nuisance estimators are used.
                            1 : T np / C np / A np
                            21: T np / C np / A p
                            22: T np / C p  / A np
                            23: T np / C p  / A p
                            3 : T p  / C np / A np
                            41: T p  / C np / A p
                            42: T p  / C p  / A np
                            43: T p  / C p  / A p")),
  make_option(c("-p", "--p"), action = "store", default = NA, type = "character",
              help = "Policy (Type B or TPB)")
)

opt = parse_args(OptionParser(option_list = option_list))

m = opt$m              # Number of clusters per simulation
r = opt$r              # Number of binary vector sampling
s = opt$s              # Scenario indicating which nuisance estimators are used
policy = opt$p         # Policy (Type B or TPB)

## Wrap TASK_ID with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

taus = 0.1*(1:5)
if(policy == "TypeB"){
  thetas = 0.3 + 0.01*0:30
  theta0 = 0.3 + 0.01*15 ## Prevent potential floating issue 0.3+0.05*15 != 0.45
  
} else if(policy == "TPB"){
  thetas <- 0.025*0:20
  theta0 = 0
  
} else {
  stop("Scenario policy (p) is incorrect. Terminate procedure.")
}

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
data = readRDS(glue("../../data/m{m}/data_id{task_ID}.rds"))


###-------- Fit estimator --------###
## Nuisance estimation scenarios
scenario_df <- tibble::tibble(
  s = c(1, 21, 22, 23, 3, 41, 42, 43),
  Tfit = c("nonpara", "nonpara", "nonpara", "nonpara", "para", "para", "para", "para"),
  Cfit = c("nonpara", "nonpara", "para",    "para",    "nonpara", "nonpara", "para", "para"),
  Afit = c("nonpara", "para",    "nonpara", "para",    "nonpara", "para",    "nonpara", "para")
)

selected <- scenario_df %>%
  dplyr::filter(s == !!s)

if (nrow(selected) == 1) {
  Tfit.method <- selected$Tfit
  Cfit.method <- selected$Cfit
  Afit.method <- selected$Afit
} else {
  stop("Scenario specified (s) is incorrect. Terminate procedure.")
}

## Covariate names
X.names = colnames(data %>% select(-c(id, Y, D, A)))

## Number of sample splints
K = 2

result = estimator(data = data,
                   X.T.names = X.names, X.C.names = X.names, X.A.names = X.names,
                   policy = policy,
                   taus = taus, thetas = thetas, theta0 = theta0, 
                   K = K, r = r, bounding = FALSE,
                   Tfit.method = Tfit.method, Cfit.method = Cfit.method, Afit.method = Afit.method)

# cat("\n")
# cat("==========================================================================\n")
# cat("==========================   Estimantion Result  =========================\n")
# cat("==========================================================================\n\n")
# 
# result$result %>%
#   mutate(across(c(est, se, PCL, PCU, UCL, UCU), function(x) round(x,5)))
# 
# cat("\n")
# 
# cat("==========================================================================\n")
# cat("==========================   Weight distribution  ========================\n")
# cat("==========================================================================\n\n")
# 
# cat("--- OR weights ---\n")
# 
# W = lapply(result$OR.w, function(l){
#   l = do.call(rbind, l)
#   melt(l, varnames = c("estimand", "theta"), value.name = "weight")
# }
# )
# 
# W = bind_rows(W, .id = "id")
# 
# W %>%
#   group_by(estimand, theta) %>%
#   summarise(
#     mean_value = mean(weight),
#     median_value = median(weight),
#     min_value = min(weight),
#     max_value = max(weight),
#     .groups = "drop"
#   ) %>% print(n=100)
# 
# 
# cat("\n--- IPCW weights ---\n")
# 
# W = lapply(result$IPCW_BC.w, function(l){
#   l = do.call(rbind, l)
#   melt(l, varnames = c("estimand", "theta"), value.name = "weight")
# }
# )
# 
# W = bind_rows(W, .id = "id")
# 
# W %>%
#   group_by(estimand, theta) %>%
#   summarise(
#     mean_value = mean(weight),
#     median_value = median(weight),
#     min_value = min(weight),
#     max_value = max(weight),
#     .groups = "drop"
#   ) %>% print(n=100)

###-------- Save simulated estimator list as Rdata --------###

## Save output
saveRDS(result$result, file = paste0("Rdata/estimate_id", task_ID,".rds"))

## Stop timer and report total run time
script.time = proc.time() - time
cat(glue("\n--- Total run time: {round(script.time[3], 1)} seconds --- \n\n"))
