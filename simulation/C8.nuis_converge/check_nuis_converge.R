time = proc.time()

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.5")
library(dplyr)
library(glue)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

library(doMC)
n.cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
n.cpus <- ifelse(is.na(n.cpus), 3, n.cpus)
registerDoMC(cores = n.cpus)
print(glue("Running in parallel with {n.cpus} cores using doMC..."))

suppressPackageStartupMessages(require(optparse))

###----------- Read arguments (simulation parameters) ---------------###
option_list = list(
  make_option(c("-m", "--m"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters per each simulation"))
)

opt = parse_args(OptionParser(option_list = option_list))

m = opt$m              # Number of clusters per simulation

## Wrap TASK_ID with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

taus = 0.1*(1:5)

###----------- Source help functions  ---------------###

code.dir = "~/research/NPSACI/code"

## Help functions for simulated data
source("../../help_simul.R")

## Help functions for nuisance functions estimation method
source(paste0(code.dir,"/help_nuis_est.R"))



###----------- Print analysis setting  ---------------###

data_train = data.sim(m)

Ni = data_train %>% 
  group_by(id) %>%
  summarise(n = n())

print("[Analysis setting]")
print(paste0("taus: ", paste(taus, collapse = ", ")))
print(paste0("m: ", nrow(Ni)))
print(paste0("N_total: ", nrow(data_train)))

print("cluster size dist'n: ")
table(Ni$n)   ## distribution of cluster size

print("observed time summary: ")
summary(data_train$Y)  ## Observed time

print("event indicator table: ")
table(data_train$D)      ## Censoring

print("treatment indicator table: ")
table(data_train$A)      ## Treatment



###-------- True nuisance functions --------###
H.true = Amodel(X = NULL, N = NULL, type="function")
F.true = Tmodel(A = NULL, X = NULL, N = NULL, type="function")
G.true = Cmodel(A = NULL, X = NULL, N = NULL, type="function")



###-------- Train nuisance estimators --------###

print("")

X.names = colnames(data_train %>% select(-c(id,Y,D,A)))

print("Parametric nuisance estimators training")

H.fit.para = H.train(data_train, X.names, "para")
F.fit.para = F.train(data_train, X.names, "para")
G.fit.para = G.train(data_train, X.names, "para")

print("Nonparametric nuisance estimators training")

H.fit.nonpara = H.train(data_train, X.names, "nonpara")
F.fit.nonpara = F.train(data_train, X.names, "nonpara")
G.fit.nonpara = G.train(data_train, X.names, "nonpara")


###-------- Evaluate estimators using new data --------###
M = 1000
data_eval = data.sim(M)
print("Evaluation using data (Parallel)")

results_list <- foreach(
  i = 1:M, 
  .combine = 'rbind', 
  .packages = c("dplyr")
) %dopar% {
    O.i = data_eval %>% filter(id == i)
    N.i = nrow(O.i)
    X.i = O.i %>% select(all_of(X.names))
    A.i = O.i$A
    
    ## True values
    H.i.true = H.true(A.i, X.i, N.i)
    F.i.true = F.true(taus, A.i, X.i, N.i)
    G.i.true = G.true(taus, A.i, X.i, N.i)
    
    ## Parametric estimate values
    H.i.hat.para = H(A.i, X.i, N.i, H.fit.para)
    F.i.hat.para = F(taus, A.i, X.i, N.i, F.fit.para)
    G.i.hat.para = G(taus, A.i, X.i, N.i, G.fit.para)
    
    ## Nonparametric estimate values
    H.i.hat.nonpara = H(A.i, X.i, N.i, H.fit.nonpara)
    F.i.hat.nonpara = F(taus, A.i, X.i, N.i, F.fit.nonpara)
    G.i.hat.nonpara = G(taus, A.i, X.i, N.i, G.fit.nonpara)
    
    ## Compute absolute differences
    h_bp = (H.i.hat.para - H.i.true)^2
    h_bn = (H.i.hat.nonpara - H.i.true)^2
    
    # Compute the column-wise vector 2-norm first (difference for each time point)
    f_bp = mean(colSums((F.i.hat.para - F.i.true)^2))
    f_bn = mean(colSums((F.i.hat.nonpara - F.i.true)^2))
    
    g_bp = mean(colSums((G.i.hat.para - G.i.true)^2))
    g_bn = mean(colSums((G.i.hat.nonpara - G.i.true)^2))
    
    # ## Compute relative differences
    # h_bp = (H.i.hat.para - H.i.true)^2 / (H.i.true + 1e-16)^2
    # h_bn = (H.i.hat.nonpara - H.i.true)^2 / (H.i.true + 1e-16)^2
    # 
    # # Compute the column-wise vector 2-norm first (difference for each time point)
    # f_bp = mean(colSums((F.i.hat.para - F.i.true)^2) / colSums(F.i.true^2 + 1e-8))
    # f_bn = mean(colSums((F.i.hat.nonpara - F.i.true)^2) / colSums(F.i.true^2 + 1e-8))
    # 
    # g_bp = mean(colSums((G.i.hat.para - G.i.true)^2) / colSums(G.i.true^2 + 1e-8))
    # g_bn = mean(colSums((G.i.hat.nonpara - G.i.true)^2) / colSums(G.i.true^2 + 1e-8))
    
    # Return a single row of the results
    data.frame(
      H_p = h_bp,
      H_np = h_bn,
      F_p = f_bp,
      F_np = f_bn,
      G_p = g_bp,
      G_np = g_bn
    )
  }

###-------- Reformat results for result.summary --------###
# This maps the parallel output back to your existing results structure
result.df = rbind(
  data.frame(bias = results_list$H_p,  method = "para",    nuis = "H"),
  data.frame(bias = results_list$H_np, method = "nonpara", nuis = "H"),
  data.frame(bias = results_list$F_p,  method = "para",    nuis = "F"),
  data.frame(bias = results_list$F_np, method = "nonpara", nuis = "F"),
  data.frame(bias = results_list$G_p,  method = "para",    nuis = "G"),
  data.frame(bias = results_list$G_np, method = "nonpara", nuis = "G")
)

result.summary = result.df %>%
  group_by(method, nuis) %>%
  summarize(bias = sqrt(mean(bias)), .groups = "drop")

## Final output
print(result.summary)

###-------- Save simulated estimator list as Rdata --------###

## Save output
saveRDS(result.summary, file = paste0("Rdata/nuis_estimate_id", task_ID,".rds"))

## Stop timer and report total run time
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))