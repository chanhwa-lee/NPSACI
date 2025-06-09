library(dplyr)
library(tidyverse)
library(glue)

## Policy specific help functions
code.dir = "~/research/NPSACI/code"
source(glue("{code.dir}/help_TypeB.R"))


###------ Cluster size (N) distribution ------###
# Define the cluster size distribution function N.dist based on ndist
N.dist <- function(m, ndist) {
  if (ndist == 1) {
    return(sample(3:5, size = m, replace = TRUE))
  }
  
  if (ndist == 2) {
    return(sample(5:20, size = m, replace = TRUE))
  }
  
  if (ndist == 3) {
    return(sample(20:50, size = m, replace = TRUE))
  }
  
  if (ndist == 4) {
    return(sample(50:100, size = m, replace = TRUE))
  }
  
  if (ndist == 5) {
    # Mimics the real data cluster size using a truncated negative binomial distribution
    # Ensures all cluster sizes are >= 2
    Nsample <- integer(0)
    while (length(Nsample) < m) {
      new <- rnbinom(ceiling(1.1 * (m - length(Nsample))), size = 1.79, mu = 19.94)
      Nsample <- c(Nsample, new[new >= 2])
    }
    return(Nsample[1:m])
  }
  
  stop("ndist must be an integer between 1 and 5")
}


###----------- Data generating process (DGP) function ---------------###

#' Data simulation and estiamnds computation based on N, X, A, T, C models
#'
#' @param N An integer. Cluster size
#' @param type A character. Indicate to generate observed data or estimand
#' @param tau A number. Time of interest
#' @param thetas A numeric vector. thetas of IPS policies. Should include 1 
#' (theta = 1 is the baseline value), otherwise 1 will be included automatically
#' @param r An integer. Number of subsampling approximation 
#' (number of random binary vector sampling)
#' @return Approximated causal estimands under the IPS policy over the theta values

## X model ##
Xmodel = function(N){
  X1 <- rnorm(N, 0, 1)
  X2 <- rnorm(N, 0, 1)
  X3 <- rnorm(N, 0, 1)
  X4 <- rnorm(N, 0, 1)
  X5 <- rnorm(N, 0, 1)
  
  X6 <- rbinom(N, 1, 0.5)
  X7 <- rbinom(N, 1, 0.5)
  X8 <- rbinom(N, 1, 0.5)
  X9 <- rbinom(N, 1, 0.5)
  X10 <- rbinom(N, 1, 0.5)
  
  Xc1 <- rnorm(1, 0, 1)
  Xc2 <- rnorm(1, 0, 1)
  Xc3 <- rnorm(1, 0, 1)
  Xc4 <- rnorm(1, 0, 1)
  Xc5 <- rnorm(1, 0, 1)
  
  X = data.frame(X1, X2, X3, X4, X5,
                 X6, X7, X8, X9, X10,
                 Xc1, Xc2, Xc3, Xc4, Xc5)
  return(X)
}

## A model ##
Amodel = function(X,N, type = "data"){
  
  sigma.b = 0.5
  fx = function(X,N) with(X, -0.1 - 0.2*X1 + 0.2*X2^2 + 0.1*I(X1>0)*X6 - 0.3*pmax(Xc1,0.5))
  
  if(type == "data"){
    b = rnorm(1,0,sigma.b)
    pi = pnorm(fx(X,N) + b)
    A <- rbinom(N, 1, pi)
    return(A)
  } 
  else if (type == "function") {
    H.true = function(a,X,N){
      f = function(b){
        pr.b <- pnorm(drop(outer(fx(X, N), b, "+")))
        hh <- (pr.b)^a * (1 - pr.b)^(1 - a)
        out <- apply(hh, 2, prod) * stats::dnorm(b, mean = 0, sd = sigma.b)
        return(out)
      }
      
      H = stats::integrate(f, lower = -Inf, upper = Inf)$value
      return(H)          
    }
    
    return(H.true)
  } 
  else if (type == "fx"){
    return(fx)
  } 
  else if (type == "sigma.b"){
    return(sigma.b)
  } 
  else {
    stop("type is incorrectly specified")
  }
  
}


## T model ##
Tmodel = function(A,X,N, type = "data"){
  
  ### Generate T from Gamma distribution - to break proportional hazard assumption
  
  shape_T <- function(A, X, N){
    g.A <- (sum(A)-A)/(N-1)
    with(X, 0.1 + 0.3*A + 0.3*sin(base::pi/2*g.A)*X1^2 + 0.1*A*g.A + 0.1*X2^2*pmax(Xc1,0.1) + 0.1*I(Xc1*Xc2 < 0.5))
  }  
  
  if(type == "data"){
    T_ <- rgamma(N, shape = shape_T(A, X, N), scale = 2)
    return(T_)
  } 
  else if (type == "function") {
    F.true <- function(taus, A, X, N){
      F.taus = sapply(taus, function(tau) pgamma(tau, shape = shape_T(A, X, N), scale = 2))
      colnames(F.taus) = taus
      return(F.taus)
    }
    return(F.true)
  } 
  else {
    stop("type is incorrectly specified")
  }
  
}

## C model ##
Cmodel = function(A,X,N, type = "data"){
  
  ### Generate C from Gamma distribution – to break proportional hazards assumption
  
  shape_C <- function(A, X, N){
    g.A <- (sum(A)-A)/(N-1)
    with(X, 0.2 + 0.5*A + 0.5*g.A*X2^2 + 0.1*pmax(X1,0.1) + 0.1*I(Xc2 < 0.5))
  } 
  
  if(type == "data"){
    C <- rgamma(N, shape = shape_C(A, X, N), scale = 2)
    return(C)
  } 
  else if (type == "function") {
    G.true <- function(taus, A, X, N){
      G.taus = sapply(taus, function(tau) 1-pgamma(tau, shape = shape_C(A, X, N), scale = 2))
      colnames(G.taus) = taus
      return(G.taus)
    }
    return(G.true)
  } 
  else {
    stop("type is incorrectly specified")
  }
  
}



### Observed data generation for one cluster ###
genCluster <- function(N){
  
  ## Step1. Covariate generation
  X = Xmodel(N)
  
  ## Step2. Treatment model
  A = Amodel(X,N)
  fx = Amodel(X,N, type = "fx")
  sigma.b = Amodel(X,N, type = "sigma.b")
  H.true = list(fx = fx, sigma.b = sigma.b)
  
  ## Step3. Survival time model
  T_ = Tmodel(A,X,N)
  C  = Cmodel(A,X,N)
  
  ## Step 4. Combine all data
  Y <- pmin(T_,C)
  D <- as.numeric(T_ <= C)
  
  data = cbind(data.frame(Y = Y, D = D, A = A), X)
  
  return(data)
  
}

### Compute estimand values for one cluster ###
computeEstimands <- function(N, policy, taus, thetas, r){
  
  ## Covariate generation
  X = Xmodel(N)
  
  ## Treatment model
  A = Amodel(X,N)
  fx = Amodel(X,N, type = "fx")
  sigma.b = Amodel(X,N, type = "sigma.b")
  H.true = list(fx = fx, sigma.b = sigma.b)
  
  ## True event time cumulative function
  F.true = Tmodel(A,X,N,type="function")
  
  ## Define w.over.H.F function, which is
  ## w/H(ai,Xi,Ni)^T F(tau|ai,Xi,Ni)
  ## where w.over.H is policy-specific function
  if (policy == "TypeB") {
    
    w.over.H.F = function(a){
      
      w.over.H.a.thetas = w.over.H.TypeB(a, X, N, thetas, H.true)  
      
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
    
  } else if (policy == "TPB") {
    
    ## Compute P.A_bar.thetas
    b.rep = rnorm(max(r,10000), 0, sigma.b)
    pr.b.rep <- pnorm(drop(outer(fx(X,N), b.rep, "+")))
    a.rep <- apply(pr.b.rep, 2, function(p) rbinom(N, 1, p))
    a.bars = colMeans(a.rep)
    P.A_bar.thetas = sapply(X = thetas, FUN = function(theta) mean(a.bars >= theta))
    names(P.A_bar.thetas) = thetas
    
    w.over.H.F = function(a){
      
      w.over.H.a.thetas = w.over.H.TPB(a, X, N, thetas, H.true, P.A_bar.thetas)
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
    
  } else {
    stop(glue("Policy specified ({policy}) is not supported."))
  }
  
  
  ## Subsampling approximation
  b.rep = rnorm(r, 0, sigma.b)
  pr.b.rep <- pnorm(drop(outer(fx(X,N), b.rep, "+")))
  a.rep <- apply(pr.b.rep, 2, function(p) rbinom(N, 1, p))
  res.thetas.taus = apply(a.rep, 2, w.over.H.F)
  
  mat.mu   = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu))
  mat.mu_1 = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu_1))
  mat.mu_0 = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu_0))
  
  df.mu = as.data.frame(mat.mu) %>%
    mutate(theta = thetas) %>%
    gather(key = "tau", value = "mu", -theta) %>%
    mutate(tau = as.numeric(tau)) %>%
    arrange(theta, tau)
  
  df.mu_1 = as.data.frame(mat.mu_1) %>%
    mutate(theta = thetas) %>%
    gather(key = "tau", value = "mu_1", -theta) %>%
    mutate(tau = as.numeric(tau)) %>%
    arrange(theta, tau)
  
  df.mu_0 = as.data.frame(mat.mu_0) %>%
    mutate(theta = thetas) %>%
    gather(key = "tau", value = "mu_0", -theta) %>%
    mutate(tau = as.numeric(tau)) %>%
    arrange(theta, tau)
  
  estimands <- df.mu %>%
    inner_join(df.mu_1, by = c("theta", "tau")) %>%
    inner_join(df.mu_0, by = c("theta", "tau")) 
  
  return(estimands)
  
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

data.sim = function(m, ndist){
  N.list = N.dist(m, ndist)
  data = lapply(N.list, genCluster)
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

estimands.sim = function(M, policy, taus, thetas, theta0, r, ndist){
  N.list = N.dist(M, ndist)
  estimands = lapply(N.list, computeEstimands, policy = policy, taus = taus, thetas = thetas, r = r)
  estimands = dplyr::bind_rows(estimands) %>% 
    group_by(theta, tau) %>% 
    summarise(mu = mean(mu),
              mu_1 = mean(mu_1),
              mu_0 = mean(mu_0),
              .groups = "drop")
  
  ### Causal effects computation using theta0 as the standard ###
  standard = estimands %>% dplyr::filter(theta == theta0)
  
  estimands = estimands %>% 
    mutate(de  = mu_1 - mu_0, 
           se_1 = mu_1 - standard$mu_1,
           se_0 = mu_0 - standard$mu_0,
           oe  = mu   - standard$mu) %>% 
    pivot_longer(cols = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe"),
                 names_to = "estimand",
                 values_to = "truth")
  
  return(estimands)
}