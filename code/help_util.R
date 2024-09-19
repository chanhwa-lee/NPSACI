library(dplyr)
library(reshape2)
library(doMC)
n.cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
n.cpus <- ifelse(is.na(n.cpus), 3, n.cpus)
registerDoMC(cores = n.cpus) 


###--- Aggregate IF function ---###
#' @input evaluated IFs over clusters {IF_i: i = 1, ..., m}
#' @intermediate
#' @output estimate and se from IFs
IF.to.est <- function(IFs, fold){
  
  d = data.frame(fold = fold, IF = IFs) %>%
    group_by(fold) %>%
    summarise(avg = mean(IF), avg2 = mean(IF^2))
  
  est = mean(d$avg)
  se = sqrt(1/length(IFs) * (mean(d$avg2) - est^2))
  
  return(c(est, se))
}


###--- Generate estimates for all estimands ---###
#' @input ifvals = list(mu, mu_1, mu_0): evaluated IFs over clusters {IF_i: i = 1, ..., m}, 
#' policy index parameters \thetas, standard index parameter \theta_0
#' @intermediate IF.to.est
#' @output estimate and se from IFs of all target estimands
estimate = function(ifvals.mu, ifvals.mu_1, ifvals.mu_0, fold, taus, thetas, theta0){
  
  ### Help function
  ifvals.to.df <- function(ifvals){
    
    ifvals.df = 
      dplyr::bind_rows(lapply(ifvals,
                       function(mat) melt(mat, varnames = c("theta", "tau")) %>% 
                         dplyr::rename(IF = value)), 
                .id = "id") %>%
      dplyr::mutate(id = as.numeric(id))
    
    return(ifvals.df)
  }
  
  df.to.est <- function(ifvals.df){
    
    m = length(fold)
    
    ifvals.df$fold = rep(fold, each = length(taus)*length(thetas))
    
    res.df = ifvals.df %>%
      dplyr::group_by(theta, tau, fold) %>%
      dplyr::summarise(avg  = mean(IF), 
                avg2 = mean(IF^2)) %>%
      dplyr::group_by(theta, tau) %>%
      dplyr::summarise(est = mean(avg),
                se  = sqrt((mean(avg2) - mean(avg)^2)/m)) %>%
      ungroup()
    
    return(res.df)
  }
  
  ### mu, mu_1, mu_0 estimates
  
  ifvals.df.mu   = ifvals.to.df(ifvals.mu)
  ifvals.df.mu_1 = ifvals.to.df(ifvals.mu_1)
  ifvals.df.mu_0 = ifvals.to.df(ifvals.mu_0)
  
  mu.est   = df.to.est(ifvals.df.mu)
  mu_1.est = df.to.est(ifvals.df.mu_1)
  mu_0.est = df.to.est(ifvals.df.mu_0)
  
  ### de   
  ifvals.df.de = 
    left_join(x = ifvals.df.mu_1, 
              y = ifvals.df.mu_0, 
              by = c("id", "theta", "tau"),
              suffix = c("_1", "_0")) %>%
    dplyr::mutate(IF = IF_1 - IF_0) %>%
    dplyr::select(-c(IF_1, IF_0))
  
  de.est = df.to.est(ifvals.df.de)
  
  ### oe
  ifvals.df.oe =
    dplyr::left_join(x = ifvals.df.mu, 
              y = ifvals.df.mu %>% filter(theta == theta0), 
              by = c("id", "tau"),
              suffix = c(".ori", ".0")) %>%
    dplyr::mutate(IF = IF.ori - IF.0) %>%
    dplyr::select(-c(IF.ori, IF.0, theta.0)) %>%
    dplyr::rename(theta = theta.ori)
  
  oe.est = df.to.est(ifvals.df.oe)
  
  ### se_1
  ifvals.df.se_1 =
    left_join(x = ifvals.df.mu_1, 
              y = ifvals.df.mu_1 %>% filter(theta == theta0), 
              by = c("id", "tau"),
              suffix = c(".ori", ".0")) %>%
    dplyr::mutate(IF = IF.ori - IF.0) %>%
    dplyr::select(-c(IF.ori, IF.0, theta.0)) %>%
    dplyr::rename(theta = theta.ori)
  
  se_1.est = df.to.est(ifvals.df.se_1)
  
  ### se_0
  ifvals.df.se_0 =
    left_join(x = ifvals.df.mu_0, 
              y = ifvals.df.mu_0 %>% filter(theta == theta0), 
              by = c("id", "tau"),
              suffix = c(".ori", ".0")) %>%
    dplyr::mutate(IF = IF.ori - IF.0) %>%
    dplyr::select(-c(IF.ori, IF.0, theta.0)) %>%
    dplyr::rename(theta = theta.ori)
  
  se_0.est = df.to.est(ifvals.df.se_0)
  
  ### te
  ifvals.df.te =
    left_join(x = ifvals.df.mu_1, 
              y = ifvals.df.mu_0 %>% filter(theta == theta0), 
              by = c("id", "tau"),
              suffix = c(".ori", ".0")) %>%
    dplyr::mutate(IF = IF.ori - IF.0) %>%
    dplyr::select(-c(IF.ori, IF.0, theta.0)) %>%
    dplyr::rename(theta = theta.ori)
  
  te.est = df.to.est(ifvals.df.te)
  
  result = list(
    mu   = mu.est,
    mu_1 = mu_1.est,
    mu_0 = mu_0.est,
    de   = de.est,
    se_1 = se_1.est,
    se_0 = se_0.est,
    oe   = oe.est,
    te   = te.est
  )
  
  result.df = bind_rows(result, .id = "estimand")
  
  return(result.df)
  
}

###--- Estimator function ---###
#' @input All data = {O.i: i=1,...,m}, 
#' time of interest \tau, if \taus = NULL, then return for observed times 
#' policy index parameter \theta and standard parameter \theta0
#' Number of splits K 
#' subsampling degree r (r = 0: no subsampling)
#' @intermediate IF, estimate
#' @output Estimator of causal estimands

estimator <- function(data, X.T.names, X.C.names, X.A.names, taus, thetas, theta0, K, r = 0){
  
  ## Record run time
  time = proc.time()
  
  ## Split dat into lists to get N for each cluster
  data.list <- split(data, f = data$id)
  
  ## setup storage
  m <- length(data.list)        # Number of clusters
  
  ifvals.mu = list()
  ifvals.mu_1 = list()
  ifvals.mu_0 = list()
  
  OR   = list()
  OR.w = list()
  
  IPCW_BC   = list()
  IPCW_BC.w = list()
  
  AUG   = list()
  AUG.w = list()
  
  fold <- sample(1:m, replace = FALSE) %% K + 1
  
  ###--- Sample Splitting Estimator ---###
  
  ## Fit treatment and outcome model at each split and evaluate IF values
  for (split in 1:K) {
    
    print(paste("   Split:", split))
    
    train.idx = which(fold != split)
    eval.idx  = which(fold == split)
    
    data_train = bind_rows(data.list[train.idx])
    
    ### Nuisance function estimation
    F.fit = F.train(data_train, X.T.names)
    G.fit = G.train(data_train, X.C.names)
    H.fit = H.train(data_train, X.A.names)
    
    
    ## First, compute OR, IPCW-BC, AUG, and their weights for eval fold clusters
    for(i in eval.idx){
      
      print(paste("     Cluster:", i))
      
      O.i = data %>% filter(id == i) %>% select(-id)
      Y.i = O.i$Y
      D.i = O.i$D
      A.i = O.i$A
      X.i = O.i %>% select(-c(Y, D, A))
      N.i = nrow(O.i)
      
      ### Compute IF value for cluster i ###
      IF.i = IF(Y.i, D.i, A.i, X.i, N.i, taus, thetas, r, F.fit, G.fit, H.fit)
      
      OR[[i]]   = IF.i$OR
      OR.w[[i]] = IF.i$OR.w
      
      IPCW_BC[[i]]   = IF.i$IPCW_BC
      IPCW_BC.w[[i]] = IF.i$IPCW_BC.w
      
      AUG[[i]]   = IF.i$AUG
      AUG.w[[i]] = IF.i$AUG.w
      
    }
    
    ## Second, compute sample-combined weights for OR, IPCW_BC, AUG in eval fold for thetas
    ## Note: expectation of all values = m_k (number of clusters in eval fold)
    OR_weight = list()
    OR_weight$mu   = Reduce("+", lapply(eval.idx, function(i) OR.w[[i]]$mu))    ## length(thetas) vector
    OR_weight$mu_1 = Reduce("+", lapply(eval.idx, function(i) OR.w[[i]]$mu_1))
    OR_weight$mu_0 = Reduce("+", lapply(eval.idx, function(i) OR.w[[i]]$mu_0))

    IPCW_BC_weight = list()
    IPCW_BC_weight$mu   = Reduce("+", lapply(eval.idx, function(i) IPCW_BC.w[[i]]$mu))
    IPCW_BC_weight$mu_1 = Reduce("+", lapply(eval.idx, function(i) IPCW_BC.w[[i]]$mu_1))
    IPCW_BC_weight$mu_0 = Reduce("+", lapply(eval.idx, function(i) IPCW_BC.w[[i]]$mu_0))

    AUG_weight = list()
    AUG_weight$mu   = Reduce("+", lapply(eval.idx, function(i) AUG.w[[i]]$mu))
    AUG_weight$mu_1 = Reduce("+", lapply(eval.idx, function(i) AUG.w[[i]]$mu_1))
    AUG_weight$mu_0 = Reduce("+", lapply(eval.idx, function(i) AUG.w[[i]]$mu_0))
    
    ## Third, re-weight OR, IPCW-BC, AUG with weights (instead of sample mean divided by m_k)
    
    m_k = length(eval.idx)

    for(i in eval.idx){

      ### New evaluated IF is the sum of re-weighted OR, IPCW-BC, AUG

      ## IF for mu at cluster i over thetas
      ## ifvals.mu  [[i]]: [length(thetas)] x[length(taus)] matrix
      ifvals.mu  [[i]] <-
        OR[[i]]$mu       +
        IPCW_BC[[i]]$mu  +
        AUG[[i]]$mu

      ## IF for mu_1 at cluster i over thetas
      ## ifvals.mu_1[[i]]: [length(thetas)] x[length(taus)] matrix
      ifvals.mu_1[[i]] <-
        OR[[i]]$mu_1      +
        IPCW_BC[[i]]$mu_1 +
        AUG[[i]]$mu_1
      
      ## IF for mu_0 at cluster i over thetas
      ## ifvals.mu_0[[i]]: [length(thetas)] x[length(taus)] matrix
      ifvals.mu_0[[i]] <-
        OR[[i]]$mu_0      +
        IPCW_BC[[i]]$mu_0 +
        AUG[[i]]$mu_0

    }

  }
  
  print("")
  print("IF for proposed nonparametric estimator computed.")
  
  ###--- Compute estimates & se from IF ---###
  
  result = estimate(ifvals.mu, ifvals.mu_1, ifvals.mu_0,
                    fold, taus, thetas, theta0)
  
  script.time = proc.time() - time
  print("")
  print(paste0("Total run time was ", script.time[3], " seconds"))
  
  return(list(result = result, 
              OR.w = OR.w, IPCW_BC.w = IPCW_BC.w))
  
}