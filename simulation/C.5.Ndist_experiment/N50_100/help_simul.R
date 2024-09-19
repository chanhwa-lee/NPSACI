library(dplyr)
library(tidyverse)

###------------------- Set rood directory ----------------------###
root.dir = "~/research/NPSACI/"

## Help functions for TypeB specific functions
source(paste0(root.dir,"/code/help_TypeB_func.R"))


###----------- Cluster size (N) distribution ---------------###

N.dist = function(m) sample(x = 50:100, size = m, replace = T)

###----------- Estimand computation function ---------------###

#' Causal estimands computation under the Incremental Propensity Score policy
#'
#' @param N An integer. Cluster size
#' @param type A character. Indicate to generate observed data or estimand
#' @param tau A number. Time of interest
#' @param thetas A numeric vector. thetas of IPS policies. Should include 1 
#' (theta = 1 is the baseline value), otherwise 1 will be included automatically
#' @param r An integer. Number of subsampling approximation 
#' (number of random binary vector sampling)
#' @return Approximated causal estimands under the IPS policy over the theta values
#' @example 
#' data = DGP(N = 5, type = "data")
#' estimand = DGP(N = 5, type = "estimands", tau = 0.9, thetas = c(0.5,1,2), r = 100)

DGP <- function(N, type = "data", taus = NULL, thetas = NULL, r = NULL){
  
  ### Observed data generation ###
  
  ## Step1. Covariate generation
  X1 <- runif(N, 0, 1)
  X2 <- runif(N, 0, 1)
  X3 <- runif(N, 0, 1)
  X4 <- runif(N, 0, 1)
  X5 <- runif(N, 0, 1)
  
  X6 <- rbinom(N, 1, 0.5)
  X7 <- rbinom(N, 1, 0.5)
  X8 <- rbinom(N, 1, 0.5)
  X9 <- rbinom(N, 1, 0.5)
  X10 <- rbinom(N, 1, 0.5)
  
  Xc1 <- runif(1, 0, 1)
  Xc2 <- runif(1, 0, 1)
  Xc3 <- runif(1, 0, 1)
  Xc4 <- runif(1, 0, 1)
  Xc5 <- runif(1, 0, 1)
  
  X = data.frame(X1, X2, X3, X4, X5,
                 X6, X7, X8, X9, X10,
                 Xc1, Xc2, Xc3, Xc4, Xc5)
  
  ## Step2. Treatment model
  sigma.b = 0.5
  b = rnorm(1,0,sigma.b)
  fx = function(X,N) with(X, 0.1 + 0.2*X1^2 + 0.2*pmax(X2,0.3)*X6 + 0.3*I(Xc1 > 0.5))
  pi = pnorm(fx(X,N) + b)
  A <- rbinom(N, 1, pi)
  
  H.true = list(fx = fx, sigma.b = sigma.b)
  
  
  ## Step3. Survival time model
  
  ### Generate T from exponential distribution
  
  shape_T <- function(A, X, N){
    g.A <- (sum(A)-A)/(N-1)
    with(X, 0.5*A + 0.4*g.A*X1 + 0.2*A*g.A + 0.2*X2 + 0.4*Xc1)
  } 
  
  T_ <- rgamma(N, shape = shape_T(A, X, N), scale = 2)
  
  F.true <- function(taus, A, X, N){
    shape = shape_T(A, X, N)
    F.taus = sapply(taus, function(tau) pgamma(tau, shape = shape, scale = 2))
    colnames(F.taus) = taus
    return(F.taus)
  }
  
  ### Generate C from exponential distribution
  
  ## Uniform (X_ij7 == 0)
  u_C <- 2*A*X1 + 0.8*X2 + 0.8*Xc1^2
  C_unif <- runif(N, 0, u_C)       
  
  ## Poisson (X_ij7 == 1)
  labmda_C <- 0.8*A*X1 + 1.0*X2 + 1.2*Xc1
  C_poi = rpois(N, labmda_C)
  
  ## C is the mixture of two survival distributions
  C <- X7 * C_unif + (1-X7) * C_poi
  
  ## Step 4. Combine all data
  Y <- pmin(T_,C)
  D <- as.numeric(T_ <= C)
  
  data = cbind(data.frame(Y = Y, D = D, A = A), X)
  
  # summary(data)
  # mean(T_ < tau)
  # mean(F_T())
  
  if(type == "data"){
    
    ### Observed data ###
    return(data)
    
  }
  
  if(type == "estimands"){
    
    #### Estimands computation over `theta` values ###
    
    ## Subsampling approximation
    b.rep = rnorm(r, 0, sigma.b)
    pr.b.rep <- pnorm(drop(outer(fx(X, N), b.rep, "+")))
    a.rep <- apply(pr.b.rep, 2, function(p) rbinom(N, 1, p))
    
    ### Define w.over.H.F function, which is
    ### w/H(ai,Xi,Ni)^T F(tau|ai,Xi,Ni)
    ### where w.over.H is policy-specific function
    
    w.over.H.F = function(a){
      
      w.over.H.a.thetas = w.over.H(a, X, N, thetas, H.true)
      F.a.taus = F.true(taus, a, X, N)
      
      ### OR function for mu(Q) over thetas
      w.over.H.F.thetas.taus   = t(w.over.H.a.thetas$mu)   %*% (F.a.taus)
      
      ### OR function for mu_1(Q) over thetas
      w.over.H.F_1.thetas.taus = t(w.over.H.a.thetas$mu_1) %*% (F.a.taus)
      
      ### OR function for mu_0(Q) over thetas
      w.over.H.F_0.thetas.taus = t(w.over.H.a.thetas$mu_0) %*% (F.a.taus)
      
      ### Result matrix
      result = list(mu   = w.over.H.F.thetas.taus, 
                    mu_1 = w.over.H.F_1.thetas.taus, 
                    mu_0 = w.over.H.F_0.thetas.taus)
      
      return(result)
      
    }
    
    # res.thetas.taus <- foreach(q = 1:r) %dopar% w.over.H.F(a.rep[,q])
    
    res.thetas.taus = apply(a.rep, 2, w.over.H.F)
    
    mat.mu   = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu))
    mat.mu_1 = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu_1))
    mat.mu_0 = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu_0))
    
    df.mu = as.data.frame(mat.mu) %>%
      rownames_to_column(var = "theta") %>%
      gather(key = "tau", value = "mu", -theta) %>%
      arrange(theta, tau)
    
    df.mu_1 = as.data.frame(mat.mu_1) %>%
      rownames_to_column(var = "theta") %>%
      gather(key = "tau", value = "mu_1", -theta) %>%
      arrange(theta, tau)
    
    df.mu_0 = as.data.frame(mat.mu_0) %>%
      rownames_to_column(var = "theta") %>%
      gather(key = "tau", value = "mu_0", -theta) %>%
      arrange(theta, tau)
    
    estimands <- df.mu %>%
      inner_join(df.mu_1, by = c("theta", "tau")) %>%
      inner_join(df.mu_0, by = c("theta", "tau")) %>%
      mutate(theta = as.numeric(theta),
             tau = as.numeric(tau))
    
    return(estimands)
    
  }
  
}

## Define true H function - used for w.over.H function
H = function(a, X, N, H.fit){
  
  fx = H.fit$fx
  sigma.b = H.fit$sigma.b
  
  f = function(b){
    pr.b <- pnorm(drop(outer(fx(X, N), b, "+")))
    hh <- (pr.b)^a * (1 - pr.b)^(1 - a)
    out <- apply(hh, 2, prod) * stats::dnorm(b, mean = 0, sd = sigma.b)
    return(out)
  }
  
  H.hat = stats::integrate(f, lower = -Inf, upper = Inf)$value
  
  return(H.hat)
}


###----------- Observed data simulation function ---------------###

#' Simulating observed data
#'
#' @param m An integer. Number of clusters
#' @return Simulated data from DGP
#' @example
#' data = data.sim(100) 
#' data = data %>% 
#'  group_by(id) %>%
#'  mutate(g.A = (sum(A) - A)/(n() - 1))

data.sim = function(m){
  N.list = N.dist(m)
  data = lapply(N.list, DGP)
  data = dplyr::bind_rows(data,.id = "id") %>% mutate(id = as.numeric(id))
  return(data)
}

###----------- Estimands simulation function ---------------###

#' Causal estimands computation under the Incremental Propensity Score policy
#'
#' @param M An integer. Number of clusters for computing estimands
#' @param tau A number. Time of interest
#' @param thetas A numeric vector. thetas of IPS policies. Should include 1 
#' (theta = 1 is the baseline value), otherwise 1 will be included automatically
#' @param r An integer. Number of subsampling approximation 
#' (number of random binary vector sampling)
#' @return Approximated causal estimands under the IPS policy over the theta values

estimands.sim = function(M, taus, thetas, theta0, r){
  N.list = N.dist(M)
  estimands = lapply(N.list, DGP, type = "estimands", taus = taus, thetas = thetas, r = r)
  estimands = dplyr::bind_rows(estimands) %>% 
    group_by(theta, tau) %>% 
    summarise(mu = mean(mu),
              mu_1 = mean(mu_1),
              mu_0 = mean(mu_0))
  
  ### Causal effects computation using theta==1 as the standard ###
  standard = estimands %>% filter(theta == theta0)
  
  estimands = estimands %>% 
    mutate(de  = mu_1 - mu_0, 
           se_1 = mu_1 - standard$mu_1,
           se_0 = mu_0 - standard$mu_0,
           oe  = mu   - standard$mu,
           te  = mu_1 - standard$mu_0) %>% 
    pivot_longer(cols = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te"),
                 names_to = "estimand",
                 values_to = "truth")
  
  return(estimands)
}
