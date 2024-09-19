time = proc.time()

###------------------- Set rood directory ----------------------###
root.dir = "~/research/NPSACI"

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(survival)
library(doMC)
library(SuperLearner)
library(randomForestSRC)
library(reshape2)
library(lme4)
library(dbarts)


suppressPackageStartupMessages(require(optparse))

###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-r", "--r"), action = "store", default = NA, type = "integer",
              help = paste0("Number of binary vector sampling for outcome reg computation")),
  make_option(c("-p", "--p"), action = "store", default = NA, type = "character",
              help = paste0("para or nonpara estimator")),
  make_option(c("-n", "--n"), action = "store", default = NA, type = "integer",
              help = paste0("maximum cluster size"))
)

opt = parse_args(OptionParser(option_list = option_list))

r = opt$r              # Number of binary vector sampling
p = opt$p              # para or nonpara estimator
n.max = opt$n          # maximum cluster size

taus = 10*(1:50)
rhos <- seq(0, 0.9, length.out = 91)
rho0 = 0

## Number of cpus (cores) for parallel computing
n.cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
n.cpus <- ifelse(is.na(n.cpus), 3, n.cpus)
registerDoMC(cores = n.cpus) 

###----------- Source help functions  ---------------###

## Help functions for estimator main functions
source(paste0(root.dir,"/code/help_util_BDD.R"))

## Help functions for TPB specific functions
source(paste0(root.dir,"/code/help_TPB_func.R"))

## Help functions for nuisance functions estimation method
if(p == "para"){
  source(paste0(root.dir,"/code/help_nuis_est_para.R"))
} else {
  source(paste0(root.dir,"/code/help_nuis_est_nonpara.R"))
}

###----------- Import dataset  ---------------###
data.ori = read.csv(paste0(root.dir,"/application/cholera_data.csv"))

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
print("[Analysis setting]")
print(paste0("n.cpus: ", n.cpus))
print(paste0("method: ", p))
print(paste0("taus: ", paste(taus, collapse = ", ")))
print(paste0("rhos: ", paste(signif(rhos, 4), collapse = ", ")))
print(paste0("r: ", r))
print(paste0("m: ", nrow(Ni)))
print(paste0("N_total: ", nrow(data)))

print("cluster size dist'n: ")
table(Ni$n)   ## distribution of cluster size

print("observed time summary: ")
summary(data$Y)  ## Observed time

print("event indicator table: ")
table(data$D)      ## Censoring

print("treatment indicator table: ")
table(data$A)      ## Treatment

###-------- Fit estimator --------###
print("")
print("Estimation start")

result = estimator(data = data, 
                   X.T.names = c("age", "dist.river"), 
                   X.C.names = c("age"), 
                   X.A.names = c("age", "dist.river"),
                   taus = taus, thetas = rhos, theta0 = rho0, K = 3, r = r)

print("Estimation result: ")

print(result$result)

###------- Weight distribution ------###
print("")
print("OR weights")

W = lapply(result$OR.w, function(l){
  l = do.call(rbind, l)
  melt(l, varnames = c("estimand", "theta"), value.name = "weight")
}
)

W = bind_rows(W, .id = "id")

W.summary = W %>%
  group_by(estimand, theta) %>%
  summarise(
    mean_value = mean(weight),
    median_value = median(weight),
    min_value = min(weight),
    max_value = max(weight)
    # Add more summary statistics as needed
  )

W.summary %>%
  filter(theta %in% c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) %>%
  print(n=100)


print("")
print("IPCW weights")

W = lapply(result$IPCW_BC.w, function(l){
  l = do.call(rbind, l)
  melt(l, varnames = c("estimand", "theta"), value.name = "weight")
}
)

W = bind_rows(W, .id = "id")

W.summary = W %>%
  group_by(estimand, theta) %>%
  summarise(
    mean_value = mean(weight),
    median_value = median(weight),
    min_value = min(weight),
    max_value = max(weight)
    # Add more summary statistics as needed
  )

W.summary %>%
  filter(theta %in% c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) %>%
  print(n=100)



###-------- Save simulated estimator list as Rdata --------###

## Wrap TASK_ID with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Save output
saveRDS(result, file = paste0("Rdata/estimate_id", task_ID,".rds"))

## Stop timer and report total run time
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))