time = proc.time()

###------------------- Set rood directory ----------------------###
root.dir = "~/research/NPSACI"

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(glue)
library(survival)
library(doMC)
library(SuperLearner)
library(randomForestSRC)
library(reshape2)
library(lme4)
library(dbarts)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(reshape2::melt)

suppressPackageStartupMessages(require(optparse))

###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-n", "--n"), action = "store", default = NA, type = "integer",
              help = paste0("maximum cluster size")),
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

n.max = opt$n          # maximum cluster size
r = opt$r              # Number of binary vector sampling
s = opt$s              # Scenario indicating which nuisance estimators are used
policy = opt$p              # Policy (Type B or TPB)

taus = 10*(1:50)

K = 3 ## Number of sample splints

## Number of cpus (cores) for parallel computing
n.cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
n.cpus <- ifelse(is.na(n.cpus), 3, n.cpus)
registerDoMC(cores = n.cpus) 

###----------- Source help functions  ---------------###

## Help functions for estimator main functions
source(paste0(root.dir,"/code/help_util.R"))

## Help functions for policy specific functions

if(policy == "TypeB"){
  source(paste0(root.dir,"/code/help_TypeB.R")) 
  thetas <- seq(0.3, 0.6, length.out = 121)
  theta0 = 0.45
  
} else if(policy == "TPB"){
  source(paste0(root.dir,"/code/help_TPB.R"))
  thetas <- seq(0, 0.5, length.out = 101)
  theta0 = 0
  
} else {
  stop("Scenario policy (p) is incorrect. Terminate procedure.")
}

## Help functions for nuisance functions estimation method
source(paste0(root.dir,"/code/help_nuis_est.R"))

###----------- Import dataset  ---------------###
data.ori = read.csv(paste0(root.dir,"/application_example/toy_data.csv"))

## Restrict data for clusters with size <= n.max for pre analysis and re-label cluster id
data = data.ori %>% 
  group_by(id) %>%
  filter(n() <= n.max, n() >= 2) %>%
  ungroup() %>%
  mutate(id = match(id, unique(id)))

Ni = data %>% 
  group_by(id) %>%
  summarise(n = n())

m = nrow(Ni)  ## Number of clusters


###----------- Print analysis setting  ---------------###
cat("\n")
cat("==========================================================================","\n")
cat("===========================   Analysis setting  ==========================","\n")
cat("==========================================================================","\n")
cat("\n")
cat(paste0("- n.cpus: ", n.cpus),"\n\n")
cat(paste0("- taus: ", paste(taus, collapse = ", ")),"\n\n")
cat(paste0("- policy: ", policy),"\n\n")
cat(paste0("- thetas: ", paste(signif(thetas, 4), collapse = ", ")),"\n\n")
cat(paste0("- r: ", r),"\n\n")
cat(paste0("- N_max: ", n.max),"\n\n")
cat(paste0("- m: ", nrow(Ni)),"\n\n")
cat(paste0("- N_total: ", nrow(data)),"\n\n")
cat(paste0("- K: ", K),"\n\n")

## Nuisance estimation scenarios
scenario_df <- tibble::tibble(
  s = c(1, 21, 22, 23, 3, 41, 42, 43),
  Tfit = c("nonpara","nonpara","nonpara","nonpara","para"   ,"para"   ,"para"   ,"para"),
  Cfit = c("nonpara","nonpara","para"   ,"para"   ,"nonpara","nonpara","para"   ,"para"),
  Afit = c("nonpara","para"   ,"nonpara","para"   ,"nonpara","para"   ,"nonpara","para")
)

selected <- scenario_df %>% dplyr::filter(s == !!s)

if (nrow(selected) == 1) {
  Tfit.method <- selected$Tfit
  Cfit.method <- selected$Cfit
  Afit.method <- selected$Afit
  
  cat("\n")
  cat(sprintf("- Nuisance estimation scenario: T %s / C %s / A %s\n",
              ifelse(Tfit.method == "para", "p", "np"),
              ifelse(Cfit.method == "para", "p", "np"),
              ifelse(Afit.method == "para", "p", "np")))
} else {
  stop("Scenario specified (s) is incorrect. Terminate procedure.")
}

cat("\n - cluster size dist'n: ")
table(Ni$n)   ## distribution of cluster size

cat("\n - observed time summary: ", "\n")
summary(data$Y)  ## Observed time

cat("\n - event indicator table: ")
table(data$D)      ## Censoring

cat("\n - event time summary: \n")
summary(data$Y[data$D == 1])  ## event time

cat("\n - censoring time summary: \n")
summary(data$Y[data$D == 0])  ## censoring time

cat("\n - treatment indicator table: ")
table(data$A)      ## Treatment



###-------- Fit estimator --------###
cat("\n")
cat("==========================================================================","\n")
cat("===========================   Estimation start  ==========================","\n")
cat("==========================================================================","\n")



result = estimator(data = data,
                   X.T.names = c("age", "dist.river"),
                   X.C.names = c("age"),
                   X.A.names = c("age", "dist.river"),
                   policy = policy,
                   taus = taus, thetas = thetas, theta0 = theta0, 
                   K = K, r = r, bounding = TRUE,
                   Tfit.method = Tfit.method, Cfit.method = Cfit.method, Afit.method = Afit.method)

cat("\n")
cat("Estimation result: ","\n")

result$result
cat("")

###------- Weight distribution ------###
cat("==========================================================================","\n")
cat("==========================   Weight distribution  ========================","\n")
cat("==========================================================================","\n")

cat("\n")
cat("OR weights","\n")

W = lapply(result$OR.w, function(l){
  l = do.call(rbind, l)
  melt(l, varnames = c("estimand", "theta"), value.name = "weight")
}
)

W = bind_rows(W, .id = "id")

W %>%
  group_by(estimand, theta) %>%
  summarise(
    mean_value = mean(weight),
    median_value = median(weight),
    min_value = min(weight),
    max_value = max(weight)
    # Add more summary statistics as needed
  ) %>% print(n=1000)


cat("\n")
cat("IPCW weights","\n")

W = lapply(result$IPCW_BC.w, function(l){
  l = do.call(rbind, l)
  melt(l, varnames = c("estimand", "theta"), value.name = "weight")
}
)

W = bind_rows(W, .id = "id")

W %>%
  group_by(estimand, theta) %>%
  summarise(
    mean_value = mean(weight),
    median_value = median(weight),
    min_value = min(weight),
    max_value = max(weight)
    # Add more summary statistics as needed
  ) %>% print(n=1000)


###-------- Save simulated estimator list as Rdata --------###

## Wrap TASK_ID with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Save output
saveRDS(result, file = paste0("Rdata/estimate_id", task_ID,".rds"))

## Stop timer and report total run time
script.time = proc.time() - time
cat(paste0("Total run time was ", script.time[3], " seconds"),"\n")