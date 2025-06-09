### Implementation of IPCW estimator from Chakladar et al. (2022) 
### 'Inverse probability weighted estimators of vaccine effects accommodating 
### partial interference and censoring'
### Code adapted from https://github.com/samrosin/interference-IPCW

time = proc.time()

###------------------- Set rood directory ----------------------###
libs <- "~/R/x86_64-pc-linux-gnu-library/4.0"
user_home_directory <- "~/research/NPSACI"

###------------------- Load libraries ----------------------###
library(optimx, lib.loc=libs)
library(msm, lib.loc=libs) #required for optimx
library(sn, lib.loc=libs) #required for optimx
library(parfm, lib.loc=libs) #parametric frailty models
library(lme4) #mixed models
library(plyr)
library(dplyr) #data manipulation
library(rootSolve, lib.loc=libs) #required for geex
library(geex, lib.loc=libs) #compute ASEs

## Wrap TASK_ID with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

###----------- Read arguments (simulation parameters) ---------------###
suppressPackageStartupMessages(require(optparse))
option_list = list(
  make_option(c("-m", "--m"), action = "store", default = NA, type = "integer",
              help = paste0("Number of clusters per each simulation"))
)

opt = parse_args(OptionParser(option_list = option_list))

m = opt$m              # Number of clusters per simulation

## Help functions for simulated data
source("../../help_simul.R")

###----------- Import dataset  ---------------###
base = readRDS(glue("../../data/m{m}/data_id{task_ID}.rds"))

Ni = base %>% 
  group_by(id) %>%
  summarise(n = n())

base <- base %>%
  dplyr::rename(group = id,
                L1_ij = X1,
                L2_ij = X2,
                L3_ij = X3,
                L4_ij = X4,
                L5_ij = X5,
                L6_ij = X6,
                L7_ij = X7,
                L8_ij = X8,
                L9_ij = X9,
                L10_ij = X10,
                Lc1_ij = Xc1,
                Lc2_ij = Xc2,
                Lc3_ij = Xc3,
                Lc4_ij = Xc4,
                Lc5_ij = Xc5,
                A_ij = A,
                Delta_ij = D,
                X_ij = Y) %>%
  group_by(group) %>%
  dplyr::mutate(k = sum(A_ij) - A_ij) %>%
  ungroup() %>%
  dplyr::mutate(sim.number = task_ID)

covariate_names = c("L1_ij", "L2_ij", "L3_ij", "L4_ij", "L5_ij", "L6_ij", "L7_ij", "L8_ij", "L9_ij", "L10_ij", "Lc1_ij", "Lc2_ij", "Lc3_ij", "Lc4_ij", "Lc5_ij")

### Time of interest
time.until = 0.4

##### specifying some simulation parameters
alphas <- c(.3,.45,.6)

#### define inv.logit function
inv.logit <- function(x){exp(x)/(1+exp(x))}

###----------- Print analysis setting  ---------------###
print("[Analysis setting]")
print(paste0("taus: ", paste(time.until, collapse = ", ")))
print(paste0("thetas: ", paste(signif(alphas, 4), collapse = ", ")))
print(paste0("m: ", m))
print(paste0("N_total: ", nrow(base)))

print("cluster size dist'n: ")
table(Ni$n)   ## distribution of cluster size

print("observed time summary: ")
summary(base$X_ij)  ## Observed time

print("event indicator table: ")
table(base$Delta_ij)      ## Censoring

print("treatment indicator table: ")
table(base$A_ij)      ## Treatment



###------------- Help functions for Chakladar IPCW estimator -------------###

##### these culminate in ase_est, a function to estimate Asymptotic Standard Error
##### using sandwich variance. relies on the two geex functions, 
##### (1) an estFUN sim.estFUN and (2) an m_estimate call inside ase_est01

sim.estFUN <- function(data, prop.model, time.until){
  d.i <- sum(1-data$Delta_ij) #number of censored obs in a group
  ind.pos <- ifelse(d.i>0, 1, 0) ###indicator of positive number of censored obs
  prop.scores <- grab_psiFUN(prop.model,data %>% select(group, A_ij, all_of(covariate_names)))
  l <- seq(from = 0, to=d.i-1) # seq for the summation
  
  #covariates in propensity model
  data.covars <- t(cbind(1, as.matrix(data %>% select(all_of(covariate_names)))))
  
  #covariates for censoring model
  l.tilde <- as.matrix(data %>% select(all_of(covariate_names)))
  
  ############# NOTE: Construction of theta ######################
  # thetas.c <- c(theta.c, theta.h, theta.r)
  # 
  # # thetas.prop <- unlist(getME(p.mixed,c('beta','theta'))) 
  # 
  # fixefs <- matrix(fixef(p.mixed),nrow=1)
  # theta.s <- (getME(p.mixed,"theta"))
  # 
  # ### Chanhwa's edit ###
  # theta.s <- ifelse(theta.s > 0, theta.s, 0.1)
  # #####################
  # 
  # thetas.prop <- c(fixefs, theta.s) 
  # 
  # ### note that grab_psiFUN wants the ranef SD, not variance, 
  # ### per the help(grab_psiFUN) page
  # muhat.0s <- mu.mat[,1]; muhat.1s <- mu.mat[,2]
  # 
  # theta <- as.vector(c(thetas.c,thetas.prop,muhat.0s, muhat.1s))
  # ## 15 (c model) + 1 (hazard) + 1 (gamma) + 16 (A model w/ intercept) + 1 (A model sigma) + 9 (alpha nums) + 9 (alpha nums)
  ################################################################
  
  idx.theta.c = (1:length(covariate_names))
  idx.theta.h = tail(idx.theta.c, n=1) + 1
  idx.theta.r = idx.theta.h + 1
  idx.fixefs  = (idx.theta.r + 1):(idx.theta.r + length(covariate_names) + 1)
  idx.theta.s = tail(idx.fixefs, n=1) + 1
  idx.muhat.0s = (idx.theta.s + 1):(idx.theta.s + length(alphas))
  idx.muhat.1s = (idx.theta.s + length(alphas) + 1):(idx.theta.s + length(alphas)*2)
  
  function(theta, alphas, time.until){
    
    l.tilde.theta <- as.vector(l.tilde %*% theta[idx.theta.c]) #fn of covs
    
    G0 <- theta[idx.theta.h]*data$X_ij ###cum hazard fn
    
    pt.4 <- ifelse(ind.pos==1,
                   sum(1/(theta[idx.theta.r]+l*theta[idx.theta.r]^2)),0) ### account for no censored obs 
    
    score.c = 
      colSums( (1-data$Delta_ij) * l.tilde) - 
      ( (d.i*theta[idx.theta.r]+1) * (colSums(G0*exp(l.tilde.theta)*l.tilde)))/
      (1+ theta[idx.theta.r]*sum(G0*exp(l.tilde.theta)))
    
    score.h <- d.i/theta[idx.theta.h]-
      (d.i*theta[idx.theta.r]+1)*sum(data$X_ij*exp(l.tilde.theta))/
      (1+ theta[idx.theta.r]*sum(G0*exp(l.tilde.theta)))
    
    score.r <- d.i/theta[idx.theta.r]-
      ((1/theta[idx.theta.r]+d.i)*sum(G0*exp(l.tilde.theta)))/
      (1+ theta[idx.theta.r]*sum(G0*exp(l.tilde.theta)))  +
      log(1+theta[idx.theta.r]*sum(G0*exp(l.tilde.theta)))/(theta[idx.theta.r]^2) -
      ind.pos*pt.4
    
    
    data$s.c <- data$Delta_ij*
      (1/(theta[idx.theta.r]*theta[idx.theta.h]*data$X_ij*exp(l.tilde.theta)+1.0))^(1/theta[idx.theta.r])
    
    fhats.gen <- function(alpha){
      ######## code to get fhat0 and fhat1
      data$alpha.denom <- alpha^data$A_ij*(1-alpha)^(1-data$A_ij)
      data$numer1 <- (1/data$alpha.denom)*data$A_ij*data$Delta_ij*(data$X_ij<time.until)
      data$numer0 <- (1/data$alpha.denom)*(1-data$A_ij)*data$Delta_ij*(data$X_ij<time.until)
      
      #propensity fn. b is random effect
      propensity.fn <- Vectorize(function(b){
        
        # fixefs <- theta[5:7]
        fixefs <- theta[idx.fixefs ]
        
        h_ij <- (t(inv.logit(fixefs %*% data.covars + b))) 
        if(anyNA(h_ij)){return(0)} #account for propensity score of 0 or 1
        if(any(h_ij>.9999999999999999)){return(0)} #account for propensity score of 1
        if(any(h_ij<1e-322)){return(0)} #account for propensity score of 0
        dn <- dnorm(x=b,mean=0,sd=theta[idx.theta.s])
        if(dn==0){return(0)} 
        return(exp(sum(
          data$A_ij*log(h_ij/alpha)+(1-data$A_ij)*log((1-h_ij)/(1-alpha))
        ))*dn)
      })
      
      #generate the propensity score 
      p.hat <- stats::integrate(f=propensity.fn,lower=-Inf,upper=Inf)$value
      
      #only compute for uncensored individuals   
      data$to.sum0 <- ifelse(data$s.c==0, 0, data$numer0/(p.hat*data$s.c))
      data$to.sum1 <- ifelse(data$s.c==0, 0, data$numer1/(p.hat*data$s.c))
      
      #wrap up
      fhat0 <- (1/nrow(data))*sum(data$to.sum0)
      fhat1 <- (1/nrow(data))*sum(data$to.sum1)
      return(c(fhat0,fhat1))
    }
    
    fhats <- lapply(X=alphas, FUN=fhats.gen)
    fhats <- matrix(unlist(fhats),ncol=2,byrow=TRUE)
    
    return(c(score.c, score.h, score.r, prop.scores(theta[c(idx.fixefs, idx.theta.s)]),
             fhats[,1] - theta[idx.muhat.0s], fhats[,2] - theta[idx.muhat.1s]))
  }
}

### Estimate the vcov using geex, using the theta values from the 
### censoring model (from parfm) and propensity model ()
### Note that the variable "group" specifies the cluster 

ase_est_call <- function(base, theta, alphas, time.until, p.mixed){
  geex_obj <- m_estimate(
    estFUN = sim.estFUN,
    data = base,
    units = "group",
    compute_roots = FALSE,
    roots = theta,
    outer_args = list(prop.model = p.mixed, time.until = time.until),
    inner_args = list(alphas = alphas, time.until = time.until),
    deriv_control = setup_deriv_control(method="simple")
  )
  
  n.theta = nrow(vcov(geex_obj))
  n.alphas = length(alphas)
  
  grab.ases <- function(a){
    if(a==0){
      ixs <- (n.theta - 2* n.alphas + 1) : (n.theta - n.alphas) 
      
    } else {
      ixs <- (n.theta - n.alphas + 1) : n.theta
    } 
    return(sqrt(diag(vcov(geex_obj)[ixs,ixs])))
  }
  
  ase_0s <- grab.ases(0); ase_1s <- grab.ases(1)
  return(matrix(c(ase_0s,ase_1s),nrow=length(alphas),ncol=2))
}



#### now estimate the ASE for both muhat.0 and muhat.1 
ase_est <- function(base, alphas, theta.c, theta.h, theta.r, p.mixed,
                    mu.mat, time.until){
  
  thetas.c <- c(theta.c, theta.h, theta.r)
  
  fixefs <- matrix(fixef(p.mixed),nrow=1)
  theta.s <- (getME(p.mixed,"theta"))
  
  theta.s <- ifelse(theta.s > 0, theta.s, 0.1)
  
  thetas.prop <- c(fixefs, theta.s) 
  
  ### note that grab_psiFUN wants the ranef SD, not variance, 
  ### per the help(grab_psiFUN) page
  muhat.0s <- mu.mat[,1]; muhat.1s <- mu.mat[,2]
  
  theta <- as.vector(c(thetas.c,thetas.prop,muhat.0s, muhat.1s))
  
  ases <- ase_est_call(base, theta = theta, alphas, time.until, p.mixed)
  
  return(ases)
}

############
#
# Component functions are:
# (1) fhat_i_aalpha implements the group-level estimator for a given value of a, alpha, and t
# (2) ipcw.eval.alpha evaluates the IPCW point estimator for a given dataset and value of alpha,
#     for both values a=0,1. It also evaluates the two ASEs by calling ase_est from ase_est.R
# (3) ipcw.eval.allalphas evaluates the IPCW estimators for a given dataset for all 
#     values of alpha. Here is where the censoring and propensity models are fit. 


##### This function returns the group-level causal estimator
fhat_i_aalpha <- function(grpdt, alpha, time.until,
                          theta.s, theta.c, p.mixed, fixefs,
                          theta.h, theta.r){
  #numerators for the estimators under a=1 and a=0
  grpdt$alpha.denom <- alpha^grpdt$A_ij*(1-alpha)^(1-grpdt$A_ij)
  grpdt$numer1 <- (1/grpdt$alpha.denom)*grpdt$A_ij*grpdt$Delta_ij*(grpdt$X_ij<time.until)
  grpdt$numer0 <- (1/grpdt$alpha.denom)*(1-grpdt$A_ij)*grpdt$Delta_ij*(grpdt$X_ij<time.until)
  
  #propensity fn. b is random effect
  propensity.fn <- Vectorize(function(b){
    
    grpdt.covars <- t(cbind(1, as.matrix(grpdt %>% select(all_of(covariate_names)))))
    
    h_ij <- (t(inv.logit(fixefs %*% grpdt.covars + b))) 
    
    if(anyNA(h_ij)){return(0)} #account for propensity score of 0 or 1
    
    if(any(h_ij>.9999999999999999)){return(0)} #account for propensity score of 1
    
    if(any(h_ij<1e-322)){return(0)} #account for propensity score of 0
    
    dn <- dnorm(x=b,mean=0,sd=sqrt(theta.s))
    
    if(dn==0){return(0)} 
    
    return(exp(sum(
      grpdt$A_ij*log(h_ij/alpha)+(1-grpdt$A_ij)*log((1-h_ij)/(1-alpha))
    ))*dn)
  })
  
  #generate the propensity score 
  p.hat <- stats::integrate(f=propensity.fn,lower=-Inf,upper=Inf)$value
  
  #### censoring model
  l.tilde <- as.matrix(grpdt %>% select(all_of(covariate_names)))
  grpdt$l.tilde.theta <- l.tilde %*% theta.c
  grpdt$s.c <- grpdt$Delta_ij*(1/(theta.r*theta.h*grpdt$X_ij*exp(grpdt$l.tilde.theta)+1.0))^(1/theta.r)
  
  #only compute for uncensored individuals   
  grpdt$to.sum0 <- ifelse(grpdt$s.c==0, 0, grpdt$numer0/(p.hat*grpdt$s.c))
  grpdt$to.sum1 <- ifelse(grpdt$s.c==0, 0, grpdt$numer1/(p.hat*grpdt$s.c))
  
  #wrap up
  fhat.0 <- (1/nrow(grpdt))*sum(grpdt$to.sum0)
  fhat.1 <- (1/nrow(grpdt))*sum(grpdt$to.sum1)
  return(c(fhat.0,fhat.1))
}

##### evaluate ipcw estimator for a specific alpha, simulated dataset called "base"
ipcw.eval.alpha <- function(base, alpha, theta.c,fixefs,theta.s,
                            theta.h, theta.r, p.mixed, time.until){
  
  base.bygroup <- split(base,base$group)
  
  grp_avg_po_est <- lapply(X=base.bygroup, FUN=fhat_i_aalpha, alpha=alpha,
                           theta.s=theta.s, theta.c=theta.c, theta.h=theta.h, theta.r=theta.r, 
                           fixefs=fixefs, time.until = time.until)
  
  grp_avg_po_est <- matrix(unlist(grp_avg_po_est), ncol=2, byrow = T)
  
  pop_avg_po_est <- matrix(apply(grp_avg_po_est, 2, mean),ncol=2)
  
  return(pop_avg_po_est)
}

#### evaluate the ipcw estimator for all levels of alpha, 1 dataset
ipcw.eval.allalphas <- function(base, alphas, time.until){
  
  pf <- parfm(Surv(X_ij,1-Delta_ij) ~ L1_ij + L2_ij + L3_ij + L4_ij + L5_ij +  
                L6_ij +  L7_ij +  L8_ij +  L9_ij +  L10_ij+  
                Lc1_ij+  Lc2_ij+  Lc3_ij+  Lc4_ij+  Lc5_ij, 
              cluster="group",
              data=base %>% select(group, X_ij, Delta_ij, all_of(covariate_names)), dist="exponential", frailty="gamma")
  
  theta.c <- as.matrix(coef(pf),ncol=1)
  theta.h <- pf[,"ESTIMATE"]["lambda"] # baseline hazard
  theta.r <- pf[,"ESTIMATE"]["theta"] # gamma frailty 
  
  print("C model fitted")
  
  p.mixed <- glmer(A_ij ~  L1_ij + L2_ij + L3_ij + L4_ij + L5_ij + 
                     L6_ij +  L7_ij +  L8_ij +  L9_ij +  L10_ij+  
                     Lc1_ij+  Lc2_ij+  Lc3_ij+  Lc4_ij+  Lc5_ij + (1|group), 
                   data=base, family=binomial)
  
  fixefs <- matrix(fixef(p.mixed),nrow=1)
  theta.s <- (getME(p.mixed,"theta"))^2 #square it b/c it returns as sd 
  theta.s <- ifelse(theta.s > 0, theta.s, 0.1)
  
  print("A model fitted")
  
  all.ipcws <- lapply(X = alphas,
                      FUN = ipcw.eval.alpha, 
                      base = base,
                      theta.c = theta.c,
                      fixefs = fixefs, 
                      theta.s = theta.s, 
                      theta.h = theta.h, 
                      theta.r = theta.r, 
                      p.mixed = p.mixed, 
                      time.until = time.until)
  
  all.ipcws <- matrix(unlist(all.ipcws),ncol=2,byrow=TRUE)
  
  print("Estimates computed")
  
  ases <- ase_est(base,alphas,theta.c,theta.h,theta.r,p.mixed,
                  mu.mat = all.ipcws, time.until)
  
  print("SEs computed")
  
  all.results <- cbind(all.ipcws, ases)
  
  colnames(all.results) <- c("mu0","mu1","ASE0","ASE1")
  
  all.results <- data.frame(all.results)
  
  all.results$theta <- alphas
  all.results$tau <- time.until
  
  mu0.results = all.results %>%
    mutate(estimand = "mu_0", est = mu0, se = ASE0) %>%
    select(estimand, theta, tau, est, se) 
  
  mu1.results = all.results %>%
    mutate(estimand = "mu_1", est = mu1, se = ASE1) %>%
    select(estimand, theta, tau, est, se) 
  
  all.results = rbind(mu1.results, mu0.results) %>%
    mutate(PCL = est - 1.96*se,
           PCU = est + 1.96*se,
           UCL = 0,
           UCU = 0)
  
  print(all.results)
  return(all.results)
}

###---------- Fit estimator ----------###

result <- ipcw.eval.allalphas(base = base, alphas = alphas, time.until=time.until)

###-------- Save simulated estimator list as Rdata --------###

## Wrap TASK_ID with as.numeric to coerce from character variable
task_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Save output
saveRDS(result, file = paste0("Rdata/estimate_id", task_ID,".rds"))

## Stop timer and report total run time
script.time = proc.time() - time
print(paste0("Total run time was ", script.time[3], " seconds"))
