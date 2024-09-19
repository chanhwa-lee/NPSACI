library(dplyr)

###--- Q (policy distribution) function ---###
#' @input Triplet (a.i, X.i, N.i), policy indices thetas
#' @intermediate possibly H or pi
#' @output Q(a.i, X.i, N.i; \theta) function for \theta \in \thetas
Q = function(a.i, X.i, N.i, thetas){
  
  ### Using outer product %o% to generate N.i x length(thetas) matrix,
  ### `pi.hat.thetas`,
  ### where each column is pi.hat.theta for theta \in thetas
  pi.hat.thetas = rep(1, N.i) %o% thetas
  colnames(pi.hat.thetas) = thetas
  
  ### Take column-product to compute Q.hat for each theta
  Q.hat.thetas = apply(a.i*pi.hat.thetas + (1-a.i)*(1-pi.hat.thetas), 2, prod)
  
  return(Q.hat.thetas)
  
}


###--- phi_Q (IF of Q) function ---###
#' @input Quadruple (A.i, X.i, N.i, a.i), policy indices thetas
#' @intermediate possibly H or pi
#' @output phi_Q(A.i, X.i, N.i; a.i; \theta) function for \theta \in \thetas
phi.Q = function(A.i, X.i, N.i, a.i, thetas){
  
  ## phi = 0 for Type B policy
  phi.thetas = matrix(0, nrow = 1, ncol = length(thetas))
  colnames(phi.thetas) = thetas
  
  return(phi.thetas)
  
}

###-----------------------------------------------------------------###
###------- Alternatives of w, phi for stability & efficiency -------###
###-----------------------------------------------------------------###

###--- w/H (weight divided by H) function ---###
#' Under conditional independence assumption
#' @input Triplet (a.i, X.i, N.i), policy indices \thetas
#' @intermediate pi.fit
#' @output w_t(a.i, X.i, N.i)/H(a.i,X.i,N.i) function for t = NULL, 0, 1
w.over.H = function(a.i, X.i, N.i, thetas, H.fit){
  
  ### Using outer product %o% to generate N.i x length(thetas) matrix,
  ### `pi.hat.thetas`,
  ### where each column is pi.hat.theta for theta \in thetas
  pi.hat.thetas = rep(1,N.i) %o% thetas        
  colnames(pi.hat.thetas) = thetas
  
  ### Take column-product to compute Q.hat for each theta
  Q.hat.thetas = apply(a.i*pi.hat.thetas + (1-a.i)*(1-pi.hat.thetas), 2, prod)
  
  ## H(a.i, X.i, N.i) estimate
  H.hat = H(a.i, X.i, N.i, H.fit)
  
  ### w/H(a,X,N) for mu(Q) over theta \in thetas
  w.over.H.thetas = 1/N.i * rep(1, N.i) %o% (Q.hat.thetas/H.hat)
  
  ### w_t/H(a,X,N) for mu_t(Q) (t = 0,1)
  ### Equation is given by
  ### [N.i^-1 I(a_ij == t)/{a_ij*pi_ij + (1-a_ij)(1-pi_ij)} 
  ###  * \prod_{l \ne j} (theta a_ij + 1- a_ij) / (theta \pi_ij + 1 - \pi_ij)]_{j=1}^{N.i}
  w.over.H_1.thetas = 
    1/N.i * 
    I(a.i == 1) / (a.i*pi.hat.thetas + (1-a.i)*(1-pi.hat.thetas)) * 
    rep(1, N.i) %o% (Q.hat.thetas/H.hat)
  
  w.over.H_0.thetas = 
    1/N.i * 
    I(a.i == 0) / (a.i*pi.hat.thetas + (1-a.i)*(1-pi.hat.thetas)) * 
    rep(1, N.i) %o% (Q.hat.thetas/H.hat)
  
  return(list(mu = w.over.H.thetas, 
              mu_1 = w.over.H_1.thetas, 
              mu_0 = w.over.H_0.thetas))
  
}



###--- phi/H (phi divided by H) function ---###
#' Under conditional independence assumption
#' @input Quadruple (A.i, X.i, N.i, a.i), policy index \theta
#' @intermediate pi.fit
#' @output phi_t(A.i, X.i, N.i; a.i)/H(a.i,X.i,N.i) function for t = NULL, 0, 1
phi.over.H = function(A.i, X.i, N.i, a.i, thetas, H.fit){
  
  ## phi = 0 for Type B policy
  phi.over.H.thetas = matrix(0, nrow = N.i, ncol = length(thetas))
  colnames(phi.over.H.thetas) = thetas
  
  return(list(mu   = phi.over.H.thetas, 
              mu_1 = phi.over.H.thetas, 
              mu_0 = phi.over.H.thetas))
  
}

###--- Outcome Regression computation (w/ Subsampling approximation) ---###
#' @input number of subsampling vectors r, 
#' cluster level data (A.i, X.i, N.i), 
#' nuisance estimators (F.fit, pi.fit)
#' policy index parameter \theta
#' @intermediate F, w.over.H, phi.over.H
#' @output approximated outcome regression part for t = -1, 0, 1 in IF using subsampling approximation (if r > 0)
#' @note Using parallel computing, need n.cpus
#' @example OutReg(A.i, X.i, N.i, taus, thetas, F.fit, pi.fit, r)

OutReg = function(A.i, X.i, N.i, taus, thetas, r = 0, F.fit, H.fit){
  
  ### Define OR.over.H function, which is
  ### {w/H(ai,Xi,Ni) + phi/H(Ai,Xi,Ni,ai)}^T F(tau|ai,Xi,Ni)
  OR.over.H.i = function(a.i){
    
    w.over.H.i.a.i.thetas = w.over.H(a.i, X.i, N.i, thetas, H.fit)
    phi.over.H.i.a.i.thetas = phi.over.H(A.i, X.i, N.i, a.i, thetas, H.fit)
    F.i.a.i.taus = F(taus, a.i, X.i, N.i, F.fit)
    
    ### OR function for mu(Q) over thetas
    OR.Weight.mu = w.over.H.i.a.i.thetas$mu   + phi.over.H.i.a.i.thetas$mu      ## [Ni] x [length(thetas)] matrix
    OR.over.H.thetas.taus   = t(OR.Weight.mu)   %*% (F.i.a.i.taus)
    
    ### OR function for mu_1(Q) over thetas
    OR.Weight.mu_1 = w.over.H.i.a.i.thetas$mu_1 + phi.over.H.i.a.i.thetas$mu_1  ## [Ni] x [length(thetas)] matrix
    OR.over.H_1.thetas.taus = t(OR.Weight.mu_1) %*% (F.i.a.i.taus)
    
    ### OR function for mu_0(Q) over thetas
    OR.Weight.mu_0 = w.over.H.i.a.i.thetas$mu_0 + phi.over.H.i.a.i.thetas$mu_0  ## [Ni] x [length(thetas)] matrix
    OR.over.H_0.thetas.taus = t(OR.Weight.mu_0) %*% (F.i.a.i.taus)
    
    ### Result matrix
    result = list(mu   = OR.over.H.thetas.taus, 
                  mu_1 = OR.over.H_1.thetas.taus, 
                  mu_0 = OR.over.H_0.thetas.taus,
                  w    = colSums(OR.Weight.mu)  ,
                  w_1  = colSums(OR.Weight.mu_1),
                  w_0  = colSums(OR.Weight.mu_0))
    
    return(result)
    
  }
  
  if(r == 0){
    
    ### if r == 0: no subsampling
    # (N.i x 2^N.i) matrix, q-th column = a.i.q
    a.i.rep = t(as.matrix(expand.grid(replicate(N.i, 0:1, simplify = FALSE))))
    
    ### Compute {w(a.i.q) + phi(a.i.q)}^T F(tau | a.i.q) for q = 1, ... , 2^N.i
    ### Over mu(Q), mu_1(Q), mu_0(Q) and over thetas and over taus
    ### Resulting matrix `res.taus.thetas` is a list with length = 2^N.i
    ### res.taus.thetas[q]: list of `mu`, `mu_1`, `mu_0`
    ### res.taus.thetas[q]$`mu`: length(taus) x length(thetas) matrix of 
    ### [{w(a.i.q; thetas[d]) + phi(a.i.q; thetas[d])}^T F(taus[t] | a.i.q)]_(t,d)
    
    res.thetas.taus <- foreach(q = 1:2^N.i) %dopar% {
      OR.over.H.i.a.i.q = OR.over.H.i(a.i.rep[,q])
      H.i.a.i.q = H(a.i.rep[,q], X.i, N.i, H.fit)
      OR.i.a.i.q = list(mu   = OR.over.H.i.a.i.q$mu   * H.i.a.i.q,
                        mu_1 = OR.over.H.i.a.i.q$mu_1 * H.i.a.i.q,
                        mu_0 = OR.over.H.i.a.i.q$mu_0 * H.i.a.i.q)
      return(OR.i.a.i.q)
    }
    
    ### Compute sum_{q} {w(a.i.q) + phi(a.i.q)}^T F(tau | a.i.q) for all a.i.q \in A(N.i)
    ### The result is length(taus) x length(thetas) matrix of 
    ### [sum_{q} {w(a.i.q; thetas[d]) + phi(a.i.q; thetas[d])}^T F(taus[t] | a.i.q)]_(t,d)
    
    OR.mu   = Reduce("+", lapply(res.thetas.taus, function(res) res$mu  ))
    OR.mu_1 = Reduce("+", lapply(res.thetas.taus, function(res) res$mu_1))
    OR.mu_0 = Reduce("+", lapply(res.thetas.taus, function(res) res$mu_0))
    
    ### For no subsampling, weights = 1
    w.thetas = rep(1,length(thetas))
    names(w.thetas) = thetas
    OR.w = OR.w_1 = OR.w_0  = w.thetas
    
    OR.all = list(mu   = OR.mu,
                  mu_1 = OR.mu_1,
                  mu_0 = OR.mu_0,
                  w   = OR.w,
                  w_1 = OR.w_1,
                  w_0 = OR.w_0)
    
    return(OR.all)
    
  }else{
    
    ### if r > 0: subsampling
    
    # r random binary vectors from A(N.i) with prob dist'n H(.)
    # a.i.q \sim H(.) <=> b.i.q \sim N(0, sigma^2), a.ij.q|b.i.q \sim Ber(expit(X.ij beta + b.i.q))
    # (N.i x r) matrix, q-th column = a.i.q
    
    a.i.rep <- H.sample(r, X.i, N.i, H.fit)
    
    ### Compute {w/H(a.i.q) + phi/H(a.i.q)}^T F(tau | a.i.q) for q = 1, ... , r
    ### Over mu(Q), mu_1(Q), mu_0(Q) and over thetas and over taus
    ### Resulting matrix `res.taus.thetas` is a list with length = r
    ### res.taus.thetas[q]: list of `mu`, `mu_1`, `mu_0`
    ### res.taus.thetas[q]$`mu`: length(taus) x length(thetas) matrix of 
    ### [{w/H(a.i.q; thetas[d]) + phi/H(a.i.q; thetas[d])}^T F(taus[t] | a.i.q)]_(t,d)
    
    res.thetas.taus <- foreach(q = 1:r) %dopar% OR.over.H.i(a.i.rep[,q])
    # If foreach (parallel computing) not working
    # res.thetas.taus <- lapply(1:r, function(q) OR.over.H.i(a.i.rep[,q]))
    
    ### Compute 1/r * sum_{q} {w(a.i.q) + phi(a.i.q)}^T F(tau | a.i.q) / f(a.i.q)
    ### where f(a.i.q) = H(a.i.q,Xi,Ni) 
    ### The result is length(taus) x length(thetas) matrix of 
    ### [1/r * sum_{q} {w/H(a.i.q; thetas[d]) + phi/H(a.i.q; thetas[d])}^T F(taus[t] | a.i.q)]_(t,d)
    
    OR.mu   = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu))
    OR.mu_1 = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu_1))
    OR.mu_0 = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu_0))
    
    ### Mean of weights used in subsampling approx, given by
    ### 1/r sum_{q} {w(a.i.q) + phi(a.i.q)}^T 1 / H(a.i.q)
    ### Note: expectation = 1
    OR.w    = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$w  ))
    OR.w_1  = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$w_1))
    OR.w_0  = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$w_0))
    
    OR.sub = list(mu   = OR.mu,
                  mu_1 = OR.mu_1,
                  mu_0 = OR.mu_0,
                  w   = OR.w,
                  w_1 = OR.w_1,
                  w_0 = OR.w_0)
    
    return(OR.sub)
    
  }
  
}

###--- IF function ---###
#' @input cluster level data O.i = (Y.i, D.i, A.i, X.i, N.i), 
#' times of interest \taus, 
#' policy indices \thetas,
#' nuisance estimators \eta = (F, G, H, w, phi),
#' subsampling degree r (r = 0: no subsampling)
#' @intermediate F, G, H, w, phi, OutReg, martingale integral
#' @output evaluated IFs over thetas by taus matrix for mu, mu_1, mu_0 

IF = function(Y.i, D.i, A.i, X.i, N.i, taus, thetas, r = 0, F.fit, G.fit, H.fit){
  
  ### Wrap up functions for convenience
  F.i.taus = F(taus, A.i, X.i, N.i, F.fit)
  
  ### G.i.Y.i = {G_ij (Y_ij)}_{j=1}^{N_i}
  G.i.Y.i = G.i_Y.i(Y.i, D.i, A.i, X.i, N.i, G.fit)
  
  ### Instead of computing H.i and w.i separately, get w.over.H
  # H.i = H(A.i, X.i, N.i, pi.fit = pi.fit)
  # w.i.thetas = w(A.i, X.i, N.i, thetas, pi.fit = pi.fit)
  w.over.H.i.thetas = w.over.H(A.i, X.i, N.i, thetas, H.fit)
  
  ### martingale integral in AUG term
  mart.i.taus = martingale(Y.i, D.i, A.i, X.i, N.i, taus, F.fit, G.fit)
  
  ### 1. Outcome Regression (OR)
  
  ## IF r >= 2**N.i, then do OR_all instead of OR_sub
  rr = ifelse(r >= 2**N.i, 0, r)
  OR.thetas.taus = OutReg(A.i, X.i, N.i, taus, thetas, rr, F.fit, H.fit)
  
  OR.i = list()
  OR.i$mu   = OR.thetas.taus$mu
  OR.i$mu_1 = OR.thetas.taus$mu_1
  OR.i$mu_0 = OR.thetas.taus$mu_0
  
  OR.w.i = list()
  OR.w.i$mu   = OR.thetas.taus$w
  OR.w.i$mu_1 = OR.thetas.taus$w_1
  OR.w.i$mu_0 = OR.thetas.taus$w_0
  
  ### 2. Inverse Probability of Censoring Weighted Bias Correction (IPCW-BC)
  IPCW_BC.i = list()
  IPCW_BC.i$mu   = t(w.over.H.i.thetas$mu)   %*% (D.i / G.i.Y.i * outer(Y.i, taus,  "<=") - F.i.taus)
  IPCW_BC.i$mu_1 = t(w.over.H.i.thetas$mu_1) %*% (D.i / G.i.Y.i * outer(Y.i, taus,  "<=") - F.i.taus)
  IPCW_BC.i$mu_0 = t(w.over.H.i.thetas$mu_0) %*% (D.i / G.i.Y.i * outer(Y.i, taus,  "<=") - F.i.taus)
  
  IPCW_BC.w.i = list()
  IPCW_BC.w.i$mu   = colSums(w.over.H.i.thetas$mu)  
  IPCW_BC.w.i$mu_1 = colSums(w.over.H.i.thetas$mu_1)
  IPCW_BC.w.i$mu_0 = colSums(w.over.H.i.thetas$mu_0)
  
  ### 3. Augmented term of censoring process martingale integral (AUG)
  AUG.i = list()
  AUG.i$mu   = t(w.over.H.i.thetas$mu)   %*% mart.i.taus
  AUG.i$mu_1 = t(w.over.H.i.thetas$mu_1) %*% mart.i.taus
  AUG.i$mu_0 = t(w.over.H.i.thetas$mu_0) %*% mart.i.taus
  
  AUG.w.i = list()
  AUG.w.i$mu   = colSums(w.over.H.i.thetas$mu)  
  AUG.w.i$mu_1 = colSums(w.over.H.i.thetas$mu_1)
  AUG.w.i$mu_0 = colSums(w.over.H.i.thetas$mu_0)
  
  return(list(OR      = OR.i,
              IPCW_BC = IPCW_BC.i,
              AUG     = AUG.i,
              OR.w      = OR.w.i,
              IPCW_BC.w = IPCW_BC.w.i,
              AUG.w     = AUG.w.i))
}
