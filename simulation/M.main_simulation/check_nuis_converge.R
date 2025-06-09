time = proc.time()

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(glue)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

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

print("Evaluation using data")

H.bias.para = numeric(M)
F.bias.para = numeric(M)
G.bias.para = numeric(M)

H.bias.nonpara = numeric(M)
F.bias.nonpara = numeric(M)
G.bias.nonpara = numeric(M)

start_time <- Sys.time()

# --- Test --- #
H.true.val = list(M)
F.true.val = list(M)
G.true.val = list(M)

H.para = list(M)
F.para = list(M)
G.para = list(M)

H.nonpara = list(M)
F.nonpara = list(M)
G.nonpara = list(M)
# ------------ #


for(i in 1:M){
  
  if (i %% (M / 10) == 0 || i == M) {
    et <- round(difftime(Sys.time(), start_time, units = "secs"), 1)
    pct <- round(100 * i / M)
    bar <- paste0("[", strrep("*", round(pct/10)), strrep("-", 10 - round(pct/10)), "]")
    print(glue("{bar} {pct}% | ET: {et}s"))
  }
  
  O.i = data_eval %>% filter(id == i)
  N.i = nrow(O.i)
  X.i = O.i %>% select(all_of(X.names))
  A.i = O.i$A
  
  ## True values
  H.true.val[[i]] = H.i.true = H.true(A.i,X.i,N.i)
  F.true.val[[i]] = F.i.true = F.true(taus, A.i,X.i,N.i)
  G.true.val[[i]] = G.i.true = G.true(taus, A.i,X.i,N.i)
  
  ## Parametric estimate values
  H.para[[i]] = H.i.hat.para = H(A.i,X.i,N.i, H.fit.para)
  F.para[[i]] = F.i.hat.para = F(taus, A.i,X.i,N.i, F.fit.para)
  G.para[[i]] = G.i.hat.para = G(taus, A.i,X.i,N.i, G.fit.para)
  
  ## Nonparametric estimate values
  H.nonpara[[i]] = H.i.hat.nonpara = H(A.i,X.i,N.i, H.fit.nonpara)
  F.nonpara[[i]] = F.i.hat.nonpara = F(taus, A.i,X.i,N.i, F.fit.nonpara)
  G.nonpara[[i]] = G.i.hat.nonpara = G(taus, A.i,X.i,N.i, G.fit.nonpara)
  
  ## Compute relative difference
  H.bias.para[i]    = abs(H.i.hat.para    - H.i.true) / (H.i.true + 1e-16)
  H.bias.nonpara[i] = abs(H.i.hat.nonpara - H.i.true) / (H.i.true + 1e-16)
  
  F.bias.para[i]    = sqrt(mean((F.i.hat.para    - F.i.true)^2)/mean(F.i.true^2 + 1e-8))
  F.bias.nonpara[i] = sqrt(mean((F.i.hat.nonpara - F.i.true)^2)/mean(F.i.true^2 + 1e-8))
  
  G.bias.para[i]    = sqrt(mean((G.i.hat.para    - G.i.true)^2)/mean(G.i.true^2 + 1e-8))
  G.bias.nonpara[i] = sqrt(mean((G.i.hat.nonpara - G.i.true)^2)/mean(G.i.true^2 + 1e-8))
}

# hist(unlist(H.true.val))
# 
# F.nonpara[[2]]
# F.true.val[[2]]
# 
# G.nonpara[[3]]
# G.true.val[[3]]


result.df = rbind(
  data.frame(bias = H.bias.para, method = "para", nuis = "H"),
  data.frame(bias = H.bias.nonpara, method = "nonpara", nuis = "H"),
  
  data.frame(bias = F.bias.para, method = "para", nuis = "F"),
  data.frame(bias = F.bias.nonpara, method = "nonpara", nuis = "F"),
  
  data.frame(bias = G.bias.para, method = "para", nuis = "G"),
  data.frame(bias = G.bias.nonpara, method = "nonpara", nuis = "G")
)

result.summary = result.df %>%
  group_by(method, nuis) %>%
  summarize(bias = mean(bias), .groups = "drop")

## Final output
print(result.summary)

ggplot(result.df, aes(x = nuis, y = bias, col = method)) +
  geom_boxplot()

###-------- Save simulated estimator list as Rdata --------###

## Save output
saveRDS(result.summary, file = paste0("Rdata/nuis_estimate_id", task_ID,".rds"))

## Stop timer and report total run time
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))