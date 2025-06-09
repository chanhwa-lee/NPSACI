library(dplyr)
library(reshape2)
library(doMC)
library(tidyr)
n.cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
n.cpus <- ifelse(is.na(n.cpus), 3, n.cpus)
registerDoMC(cores = n.cpus) 


###--- Estimator function ---###
#' @input All data = {O.i: i=1,...,m}, 
#' time of interest \tau, if \taus = NULL, then return for observed times 
#' policy index parameter \theta and standard parameter \theta0
#' Number of splits K 
#' subsampling degree r (r = 0: no subsampling)
#' @intermediate IF, estimate
#' @output Estimator of causal estimands

estimator <- function(data, 
                      X.T.names, X.C.names, X.A.names, 
                      policy,
                      taus, thetas, theta0, 
                      K = 2, r = 100, bounding = TRUE,
                      Tfit.method = "nonpara", Cfit.method = "nonpara", Afit.method = "nonpara"){
  
  ## Record run time
  start_time = Sys.time()
  
  ## Number of clusters
  m <- dplyr::n_distinct(data %>% pull(id))        
  
  cat("\n==========================================================================\n")
  cat("===========================   Analysis setting  ==========================\n")
  cat("==========================================================================\n\n")
  
  cat(paste0("- taus: ", paste(taus, collapse = ", ")),"\n\n")
  cat(paste0("- policy: ", policy),"\n\n")
  cat(paste0("- thetas: ", paste(signif(thetas, 4), collapse = ", ")),"\n\n")
  cat(paste0("- r: ", r),"\n\n")
  cat(paste0("- K: ", K),"\n\n")
  cat(glue("- Nuisance estimation method: T {Tfit.method} / C {Cfit.method} / A {Afit.method}"),"\n\n")
  
  cat(paste0("- m: ", m),"\n\n")
  cat(paste0("- N_total: ", nrow(data)),"\n\n")
  
  cat("- cluster size dist'n: ")
  print(table(data %>% group_by(id) %>% summarise(n = n()) %>% pull(n)))
  
  cat("\n - observed time summary: \n")
  print(summary(data$Y))
  
  cat("\n - event indicator table: ")
  print(table(data$D))
  
  cat("\n - event time summary: \n")
  print(summary(data$Y[data$D == 1]))
  
  cat("\n - censoring time summary: \n")
  print(summary(data$Y[data$D == 0]))
  
  cat("\n - treatment indicator table: ")
  print(table(data$A))
  
  
  cat("\n==========================================================================\n")
  cat("===========================   Estimation start  ==========================\n")
  cat("==========================================================================\n")
  
  ## setup storage for influence functions
  ifvals.mu = ifvals.mu_1 = ifvals.mu_0 = list()
  
  ## setup storage for parts of IF
  OR = IPCW_BC = AUG = list()
  
  ## setup storage for IPCweights (for bounding)
  OR.w = IPCW_BC.w = AUG.w = list()
  
  ###--- Cross-fitting Estimator ---###
  
  ## Cross-fitting folds
  fold <- sample(1:m, replace = FALSE) %% K + 1
  
  ## Fit treatment and outcome model at each split and evaluate IF values
  for (split in 1:K) {
    
    cat(glue("\n\n --- Split:{split} ---\n\n"))
    
    train.ids = which(fold != split)
    eval.ids  = which(fold == split)
    
    data_train = data %>% filter(id %in% train.ids)
    
    ## Nuisance function estimation
    F.fit = F.train(data_train, X.T.names, Tfit.method)
    G.fit = G.train(data_train, X.C.names, Cfit.method)
    H.fit = H.train(data_train, X.A.names, Afit.method)
    
    ## Number of clusters in fold k
    m_k = length(eval.ids)
    
    ## Evaluation start time at fold k
    eval_start_time = Sys.time()
    progress_marks <- floor(seq(0.1, 1, by = 0.1) * m_k)
    mark_idx <- 1
    
    ## First, compute OR, IPCW-BC, AUG, and their weights for eval fold clusters
    cat("\n")
    for(idx in 1:m_k){
      
      i = eval.ids[idx]
      
      ## Show progress bar for every 10%
      if (mark_idx <= length(progress_marks) && idx == progress_marks[mark_idx]) {
        pct <- mark_idx * 10
        et <- round(difftime(Sys.time(), eval_start_time, units = "secs"), 1)
        bar <- paste0("[", strrep("*", mark_idx), strrep("-", 10 - mark_idx), "]")
        cat(glue("   {bar} {pct}% | ET: {et}s"), "\n")
        mark_idx <- mark_idx + 1
      }
      
      O.i = data %>% filter(id == i) %>% select(-id)
      Y.i = O.i$Y
      D.i = O.i$D
      A.i = O.i$A
      X.i = O.i %>% select(-c(Y, D, A))
      N.i = nrow(O.i)
      
      ## Compute IF value for cluster i
      if (policy == "TypeB") {
        IF.i = IF.TypeB(Y.i, D.i, A.i, X.i, N.i, taus, thetas, r, F.fit, G.fit, H.fit)  
        
      } else if (policy == "TPB"){
        IF.i = IF.TPB(Y.i, D.i, A.i, X.i, N.i, taus, thetas, r, F.fit, G.fit, H.fit)  
        
      } else {
        stop(glue("Policy specified ({policy}) is not supported."))
      }
      
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
    OR_weight$mu   = Reduce("+", lapply(eval.ids, function(i) OR.w[[i]]$mu))    ## length(thetas) vector
    OR_weight$mu_1 = Reduce("+", lapply(eval.ids, function(i) OR.w[[i]]$mu_1))
    OR_weight$mu_0 = Reduce("+", lapply(eval.ids, function(i) OR.w[[i]]$mu_0))
    
    IPCW_BC_weight = list()
    IPCW_BC_weight$mu   = Reduce("+", lapply(eval.ids, function(i) IPCW_BC.w[[i]]$mu))
    IPCW_BC_weight$mu_1 = Reduce("+", lapply(eval.ids, function(i) IPCW_BC.w[[i]]$mu_1))
    IPCW_BC_weight$mu_0 = Reduce("+", lapply(eval.ids, function(i) IPCW_BC.w[[i]]$mu_0))
    
    AUG_weight = list()
    AUG_weight$mu   = Reduce("+", lapply(eval.ids, function(i) AUG.w[[i]]$mu))
    AUG_weight$mu_1 = Reduce("+", lapply(eval.ids, function(i) AUG.w[[i]]$mu_1))
    AUG_weight$mu_0 = Reduce("+", lapply(eval.ids, function(i) AUG.w[[i]]$mu_0))
    
    ## Third, re-weight OR, IPCW-BC, AUG with weights (instead of sample mean divided by m_k)
    for(i in eval.ids){
      
      ## New evaluated IF is the sum of re-weighted OR, IPCW-BC, AUG
      
      ## IF for mu at cluster i over thetas
      ## ifvals.mu  [[i]]: [length(thetas)] x[length(taus)] matrix
      ifvals.mu[[i]] <-
        OR[[i]]$mu +
        (I(bounding) * m_k / IPCW_BC_weight$mu + I(!bounding)) * IPCW_BC[[i]]$mu +
        (I(bounding) * m_k / AUG_weight$mu     + I(!bounding)) * AUG[[i]]$mu
      
      ## IF for mu_1 at cluster i over thetas
      ## ifvals.mu_1[[i]]: [length(thetas)] x[length(taus)] matrix
      ifvals.mu_1[[i]] <-
        OR[[i]]$mu_1 +
        (I(bounding) * m_k / IPCW_BC_weight$mu_1 + I(!bounding)) * IPCW_BC[[i]]$mu_1 +
        (I(bounding) * m_k / AUG_weight$mu_1     + I(!bounding)) * AUG[[i]]$mu_1
      
      ## IF for mu_0 at cluster i over thetas
      ## ifvals.mu_0[[i]]: [length(thetas)] x[length(taus)] matrix
      ifvals.mu_0[[i]] <-
        OR[[i]]$mu_0 +
        (I(bounding) * m_k / IPCW_BC_weight$mu_0 + I(!bounding)) * IPCW_BC[[i]]$mu_0 +
        (I(bounding) * m_k / AUG_weight$mu_0     + I(!bounding)) * AUG[[i]]$mu_0
      
    }
    
  }
  
  cat(glue("\n\n--- IF for proposed nonparametric estimator computed --- ",
           "ET: {round(difftime(Sys.time(), start_time, units = 'secs'), 1)} seconds\n\n"))
  
  
  ###--- Compute estimate & SE & UCB from IF ---###
  
  ## Help function: Change list of IFs to long format data frame
  ifvals.to.df <- function(ifvals){
    
    ifvals.df = 
      dplyr::bind_rows(lapply(ifvals,
                              function(mat) melt(mat, varnames = c("theta", "tau")) %>% 
                                dplyr::rename(IF = value)), 
                       .id = "id") %>%
      dplyr::mutate(id = as.numeric(id))
    
    return(ifvals.df)
  }
  
  ifvals.df = list()
  
  ifvals.df[["mu"]]   = ifvals.to.df(ifvals.mu)
  ifvals.df[["mu_1"]] = ifvals.to.df(ifvals.mu_1)
  ifvals.df[["mu_0"]] = ifvals.to.df(ifvals.mu_0)
  
  
  ifvals.df[["de"]] = 
    left_join(x = ifvals.df[["mu_1"]], 
              y = ifvals.df[["mu_0"]], 
              by = c("id", "theta", "tau"),
              suffix = c("_1", "_0")) %>%
    dplyr::mutate(IF = IF_1 - IF_0) %>%
    dplyr::select(-c(IF_1, IF_0))
  
  ifvals.df[["oe"]] =
    dplyr::left_join(x = ifvals.df[["mu"]], 
                     y = ifvals.df[["mu"]] %>% filter(near(theta,theta0)), 
                     by = c("id", "tau"),
                     suffix = c(".ori", ".0")) %>%
    dplyr::mutate(IF = IF.ori - IF.0) %>%
    dplyr::select(-c(IF.ori, IF.0, theta.0)) %>%
    dplyr::rename(theta = theta.ori)
  
  ifvals.df[["se_1"]] =
    left_join(x = ifvals.df[["mu_1"]], 
              y = ifvals.df[["mu_1"]] %>% filter(near(theta,theta0)), 
              by = c("id", "tau"),
              suffix = c(".ori", ".0")) %>%
    dplyr::mutate(IF = IF.ori - IF.0) %>%
    dplyr::select(-c(IF.ori, IF.0, theta.0)) %>%
    dplyr::rename(theta = theta.ori)
  
  ifvals.df[["se_0"]] =
    left_join(x = ifvals.df[["mu_0"]], 
              y = ifvals.df[["mu_0"]] %>% filter(near(theta,theta0)), 
              by = c("id", "tau"),
              suffix = c(".ori", ".0")) %>%
    dplyr::mutate(IF = IF.ori - IF.0) %>%
    dplyr::select(-c(IF.ori, IF.0, theta.0)) %>%
    dplyr::rename(theta = theta.ori)
  
  ## Help function: Compute estimate, SE, Uniform Confidence Bound using multiplier bootstrap
  ifvals.to.estimate = function(ifvals.df){
    
    df.estimate = data.frame()
    
    for(tau.idx in 1:length(taus)){
      
      tau = taus[tau.idx]
      
      ## IF matrix at tau: (n.thetas) x m
      ifvals.tau = ifvals.df %>%
        dplyr::filter(near(tau,taus[tau.idx])) %>%     ## To avoid floating point issue
        dplyr::select(id,theta,IF) %>%
        pivot_wider(names_from = id, values_from = IF) %>%
        dplyr::select(-theta) %>%
        as.matrix()
      
      rownames(ifvals.tau) = thetas
      
      ## estimate replicated mat: (n.thetas) x m
      est.tau = rowMeans(ifvals.tau)
      est.mat.tau = matrix(rep(est.tau, m), ncol = m, byrow = FALSE)
      rownames(est.mat.tau) = thetas
      
      ## se replicated mat: (n.thetas) x m
      se.tau = apply(ifvals.tau, 1, function(ifval) sd(ifval, na.rm = TRUE) / sqrt(m-1))
      se.mat.tau = matrix(rep(se.tau, m), ncol = m, byrow = FALSE)
      rownames(se.mat.tau) = thetas
      
      ## Standarized IF at tau: (n.thetas) x m
      ifvals.tau.tilde = ifelse(se.mat.tau != 0, (ifvals.tau - est.mat.tau) / se.mat.tau, 0) ## To avoid division by 0
      
      ## Rademacher Multiplier Boostrap process for Uniform Confidence Bound 
      n.mult.boot = 10000
      mult.mat = matrix(2*rbinom(m*n.mult.boot, 1, 0.5)-1, ncol = m)
      
      ## sup_\theta \sum_i \xi_i^(b)* \tilde{\hat{\varphi}} / m, b=1,...,n.mult.boot
      max.mult.boot.process.tau = apply(ifvals.tau.tilde %*% t(mult.mat) / m, 2, function(col) max(abs(col)))
      
      ## Take 95% quantile of Multiplier Boostrap process values to get Calpha (critical value)
      calpha.tau = quantile(max.mult.boot.process.tau, 0.95, na.rm = TRUE)
      
      ## Final output
      df.estimate.tau = data.frame(
        theta = thetas,
        tau = tau,
        est = est.tau,
        se = se.tau,
        PCL = est.tau - 1.96*se.tau,         ## `P`oint-wise Wald 95% `C`I `L`ower limit
        PCU = est.tau + 1.96*se.tau,         ## `P`oint-wise Wald 95% `C`I `U`pper limit
        UCL = est.tau - calpha.tau * se.tau, ## `U`niform 95% `C`onfidence band `L`ower limit
        UCU = est.tau + calpha.tau * se.tau  ## `U`niform 95% `C`onfidence band `U`pper limit
      )
      rownames(df.estimate.tau) = NULL
      
      df.estimate = rbind(df.estimate, df.estimate.tau)
      
    }
    
    return(df.estimate)
  }
  
  result = list()
  
  for(estimand in c("mu", "mu_1", "mu_0", "de", "oe", "se_1", "se_0")){
  
    est_start_time = Sys.time()
    result[[estimand]] = ifvals.to.estimate(ifvals.df[[estimand]])
    cat(glue("\n--- Estimator computed: {estimand} --- ",
             "ET: {round(difftime(Sys.time(), est_start_time, units = 'secs'), 1)} seconds --- \n\n"))
    
  }

  ## Combine all into one dataframe
  result = dplyr::bind_rows(result, .id = "estimand")
  
  cat(glue("\n--- Total run time: {round(difftime(Sys.time(), start_time, units = 'secs'), 1)} seconds --- \n\n"))
  
  return(list(result = result, 
              OR.w = OR.w, IPCW_BC.w = IPCW_BC.w))
  
}