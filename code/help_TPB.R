library(dplyr)

###--- w/H (weight divided by H) function ---###
#' Under conditional independence assumption
#' @input Triplet (a.i, X.i, N.i), policy indices \thetas
#' @intermediate pi.fit
#' @output w_t(a.i, X.i, N.i)/H(a.i,X.i,N.i) function for t = NULL, 0, 1

w.over.H.TPB = function(a.i, X.i, N.i, thetas, H.fit, P.A_bar.thetas){
  
  ### Q/H(a.i, X.i, N.i; theta) = I(a.i_bar >= theta) / P(A_bar >= theta)
  Q.over.H.thetas = I(mean(a.i) >= thetas) / P.A_bar.thetas
  
  ### w/H(a,X,N) for mu(Q) over theta \in thetas
  w.over.H.thetas = 1/N.i * rep(1, N.i) %o% Q.over.H.thetas
  
  ### w_t/H(a,X,N) for mu_t(Q) (t = 0,1)
  ### Equation is given by
  ### [N.i^-1 I(a_ij == t) Q(a_i(-j))/H(a_i)]_{j=1}^{N.i}
  ### Here, 
  ### Q(a_i(-j))/H(a_i) = {Q(a_i) + Q(1-a_ij, a_i(-j))} / H(a_i)
  ### = Q/H(a_i) + Q(1-a_ij,a_i(-j))/H(a_i)
  ### Thus,
  ### w_t/H(a,X,N) = [
  ###   N.i^-1 I(a_ij == t) *
  ###   {Q/H(a_i) + Q(1-a_ij,a_i(-j)) / H(a_i)}
  ### ]_{j=1}^{N.i}
  
  ### Compute [Q(1-a_ij,a_i(-j))]_{j=1}^{N.i} where
  ### Q(a.i, X.i, N.i; theta) = I(a.i_bar >= theta) / P(A_bar >= theta)
  ### and divide it by H(a.i) to get 
  ### Q.over.H.thetas.star: N.i x length(thetas) matrix
  Q.thetas.star = sapply(1:N.i, function(j){
    a.i[j] = 1 - a.i[j]
    Q.a.ij.star = I(mean(a.i) >= thetas) * H(a.i, X.i, N.i, H.fit) / P.A_bar.thetas
    return(Q.a.ij.star)
  })
  
  H.a.i = H(a.i, X.i, N.i, H.fit)
  
  Q.over.H.thetas.star = t(matrix(Q.thetas.star, nrow = length(thetas))) / H.a.i
  
  colnames(Q.over.H.thetas.star) = thetas
  
  w.over.H_1.thetas = 
    1/N.i * I(a.i == 1) *
    (rep(1, N.i) %o% Q.over.H.thetas + Q.over.H.thetas.star)
  
  w.over.H_0.thetas = 
    1/N.i * I(a.i == 0) *
    (rep(1, N.i) %o% Q.over.H.thetas + Q.over.H.thetas.star)
  
  return(list(mu   = w.over.H.thetas, 
              mu_1 = w.over.H_1.thetas, 
              mu_0 = w.over.H_0.thetas))
}


###--- phi/H (phi divided by H) function ---###
#' Under conditional independence assumption
#' @input Quadruple (A.i, X.i, N.i, a.i), policy index \theta
#' @intermediate pi.fit
#' @output phi_t(A.i, X.i, N.i; a.i)/H(a.i,X.i,N.i) function for t = NULL, 0, 1
phi2.over.H = function(A.i, X.i, N.i, a.i, thetas, H.fit, P.A_bar.thetas){
  
  ### phiQ2/H(A.i, X.i, N.i; a.i ; theta) 
  ### = - I(a.i_bar >= theta) * I(A.i_bar>=a.i_bar) / P(A_bar >= theta)^2
  phiQ2.over.H.thetas = - I(mean(a.i) >= thetas) * I(mean(A.i) >= thetas) / P.A_bar.thetas^2
  
  ### phi/H(a,X,N) for mu(Q) over theta \in thetas
  phi2.over.H.thetas = 1/N.i * rep(1, N.i) %o% phiQ2.over.H.thetas
  
  ### phi_t/H(a,X,N) for mu_t(Q) (t = 0,1)
  ### Equation is given by
  ### [N.i^-1 I(a_ij == t) phiQ(a_i(-j))/H(a_i)]_{j=1}^{N.i}
  ### Here, 
  ### phiQ(a_i(-j))/H(a_i) = {phiQ(a_i) + phiQ(1-a_ij, a_i(-j))} / H(a_i)
  ### = phiQ/H(a_i) + phiQ(1-a_ij,a_i(-j)) / H(a_i).
  ### Thus,
  ### phi_t/H(a,X,N) = [
  ###   N.i^-1 I(a_ij == t) *
  ###   {phiQ/H(a_i) + phiQ(1-a_ij,a_i(-j)) / H(a_i)
  ### ]_{j=1}^{N.i}
  
  ### Compute [phiQ(1-a_ij,a_i(-j))]_{j=1}^{N.i} where
  ### phiQ2(a.i, X.i, N.i; theta) 
  ### = - I(a.i_bar >= theta) * I(A.i_bar>=a.i_bar) * H(a.i, X.i, N.i) / P(A_bar >= theta)^2
  phiQ2.thetas.star = sapply(1:N.i, function(j){
    a.i[j] = 1 - a.i[j]
    phiQ2.a.ij.star = - I(mean(a.i) >= thetas)*I(mean(A.i) >= thetas) * H(a.i, X.i, N.i, H.fit) / P.A_bar.thetas^2
    return(phiQ2.a.ij.star)
  })
  
  H.a.i = H(a.i, X.i, N.i, H.fit)
  
  phiQ2.over.H.thetas.star = t(matrix(phiQ2.thetas.star, nrow = length(thetas))) / H.a.i
  
  colnames(phiQ2.over.H.thetas.star) = thetas
  
  phi2.over.H_1.thetas = 
    1/N.i * I(a.i == 1) *
    (rep(1, N.i) %o% phiQ2.over.H.thetas + phiQ2.over.H.thetas.star)
  
  phi2.over.H_0.thetas = 
    1/N.i * I(a.i == 0) *
    (rep(1, N.i) %o% phiQ2.over.H.thetas + phiQ2.over.H.thetas.star)
  
  return(list(mu =   phi2.over.H.thetas, 
              mu_1 = phi2.over.H_1.thetas, 
              mu_0 = phi2.over.H_0.thetas))
  
}

### It is noted that when implying Subsampling approximation to compute Outcome Regression,
### the usual subsampling would not work well since 
### phi_Q(A.i, X.i, N.i; a.i ; theta) 
### = I(a.i_bar >= theta) * {I(A.i=a.i)P(A_bar >= theta)-I(A.i_bar>=a.i_bar)H(a.i, X.i, N.i)} / P(A_bar >= theta)^2
### is in different scale (much larger) when a.i == A.i compared to other a.i's.
### i.e, phi_Q is positive if a.i == A.i and negative otherwise.

### Let phi_Q = phi_Q1 + phi_Q2
### phi_Q1(A.i, X.i, N.i; a.i ; theta) 
### = I(a.i_bar >= theta) * I(A.i=a.i) / P(A_bar >= theta)
### phi_Q2(A.i, X.i, N.i; a.i ; theta) 
### = - I(a.i_bar >= theta) * I(A.i_bar>=a.i_bar) * H(a.i, X.i, N.i) / P(A_bar >= theta)^2
### Note that phi_Q1 is nonzero only if a.i == A.i
### Then,
### OR = sum_{a.i \in A(N.i)} {w(a.i) + phi2(a.i)}^T F(tau | a.i) + sum_{a.i} phi1(a.i)^T F(tau | a.i)
###    = // + phi1(A.i)^T F(tau | A.i)
### Subsampling approximation only applies to sum_{a.i \in A(N.i)} {w(a.i) + phi2(a.i)}^T F(tau | a.i)

### Therefore, subsampling approximation here is
### OR = sum_{a.i \in A(N.i)} {w(a.i) + phi(a.i)}^T F(tau | a.i)
###    =\approx= phi1(A.i)^T F(tau | A.i) + 1/r * sum_{q} {w/H(a.i.q) + phi1/H(a.i.q)}^T F(tau | a.i.q)
### where a.i.q \sim H(.) 


###--- Outcome Regression computation (w/ Subsampling approximation) ---###
#' @input number of subsampling vectors r, 
#' cluster level data (A.i, X.i, N.i), 
#' nuisance estimators (F.fit, pi.fit)
#' policy index parameter \theta
#' @intermediate F, w.over.H, phi.over.H
#' @output approximated outcome regression part for t = -1, 0, 1 in IF using subsampling approximation (if r > 0)
#' @note Using parallel computing, need n.cpus
#' @example OutReg(A.i, X.i, N.i, taus, thetas, F.fit, pi.fit, r)

OutReg.TPB = function(A.i, X.i, N.i, taus, thetas, r = 0, F.fit, H.fit, P.A_bar.thetas){
  
  ##------ I. Compute part 1 ------##
  ### First, compute summation wrt phiQ1 which is not approximated by subsampling
  ### phi_Q1(A.i, X.i, N.i; A.i ; theta) = I(A.i_bar >= theta) / P(A_bar >= theta)
  phi.Q.1.i.thetas = I(mean(A.i) >= thetas) / P.A_bar.thetas
  F.i.taus = F(taus, A.i, X.i, N.i, F.fit)
  
  ### For mu(Q),
  ### [Ni^-1 phi_Q1(Ai, Xi, Ni; Ai; theta)]_{j=1}^{Ni}^T F(tau|Ai, Xi, Ni) 
  OR.mu.part1 = t(1/N.i * rep(1, N.i) %o% phi.Q.1.i.thetas) %*% F.i.taus
  
  ### For mu_t(Q),
  ### [Ni^-1 I(A_ij == t) phi_Q1(Ai, Xi, Ni; Ai; theta)]_{j}^T F(tau|Ai, Xi, Ni)
  ### + 
  ### [Ni^-1 I(1-A_ij == t) phi_Q1(Ai, Xi, Ni; Ai; theta)]_{j}^T [F_j(tau|1-A_ij, Ai(-j), Xi, Ni)]_{j}
  
  F.i.taus.star = sapply(1:N.i, function(j){
    A.i[j] = 1-A.i[j]
    F(taus, A.i, X.i, N.i, F.fit)[j,]
  })
  
  F.i.taus.star = t(matrix(F.i.taus.star, nrow = length(taus)))
  
  colnames(F.i.taus.star) = taus
  
  OR.mu_1.part1 = 
    t(1/N.i * I(A.i == 1) *(rep(1, N.i) %o% phi.Q.1.i.thetas)) %*% F.i.taus +
    t(1/N.i * I(1-A.i == 1) *(rep(1, N.i) %o% phi.Q.1.i.thetas)) %*% F.i.taus.star
  
  OR.mu_0.part1 = 
    t(1/N.i * I(A.i == 0) *(rep(1, N.i) %o% phi.Q.1.i.thetas)) %*% F.i.taus +
    t(1/N.i * I(1-A.i == 0) *(rep(1, N.i) %o% phi.Q.1.i.thetas)) %*% F.i.taus.star
  
  
  ##------ II. Compute part 2 ------##
  
  ### Define OR.over.H function, which is
  ### {w/H(ai,Xi,Ni) + phi2/H(Ai,Xi,Ni,ai)}^T F(tau|ai,Xi,Ni)
  OR.over.H.i = function(a.i){
    
    w.over.H.i.a.i.thetas = w.over.H.TPB(a.i, X.i, N.i, thetas, H.fit, P.A_bar.thetas)
    phi2.over.H.i.a.i.thetas = phi2.over.H(A.i, X.i, N.i, a.i, thetas, H.fit, P.A_bar.thetas)
    F.i.a.i.taus = F(taus, a.i, X.i, N.i, F.fit)
    
    ### OR function for mu(Q) over thetas
    OR.Weight.mu = w.over.H.i.a.i.thetas$mu   + phi2.over.H.i.a.i.thetas$mu      ## [Ni] x [length(thetas)] matrix
    OR.over.H.thetas.taus   = t(OR.Weight.mu)   %*% (F.i.a.i.taus)
    
    ### OR function for mu_1(Q) over thetas
    OR.Weight.mu_1 = w.over.H.i.a.i.thetas$mu_1 + phi2.over.H.i.a.i.thetas$mu_1  ## [Ni] x [length(thetas)] matrix
    OR.over.H_1.thetas.taus = t(OR.Weight.mu_1) %*% (F.i.a.i.taus)
    
    ### OR function for mu_0(Q) over thetas
    OR.Weight.mu_0 = w.over.H.i.a.i.thetas$mu_0 + phi2.over.H.i.a.i.thetas$mu_0  ## [Ni] x [length(thetas)] matrix
    OR.over.H_0.thetas.taus = t(OR.Weight.mu_0) %*% (F.i.a.i.taus)
    
    ### Result matrix
    result = list(mu   = OR.over.H.thetas.taus, 
                  mu_1 = OR.over.H_1.thetas.taus, 
                  mu_0 = OR.over.H_0.thetas.taus,
                  w    = phi.Q.1.i.thetas + colSums(OR.Weight.mu)  ,
                  w_1  = phi.Q.1.i.thetas + colSums(OR.Weight.mu_1),
                  w_0  = phi.Q.1.i.thetas + colSums(OR.Weight.mu_0))
    
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
    
    OR.mu.part2   = Reduce("+", lapply(res.thetas.taus, function(res) res$mu  ))
    OR.mu_1.part2 = Reduce("+", lapply(res.thetas.taus, function(res) res$mu_1))
    OR.mu_0.part2 = Reduce("+", lapply(res.thetas.taus, function(res) res$mu_0))
    
    ### For no subsampling, weights = 1
    w.thetas = rep(1,length(thetas))
    names(w.thetas) = thetas
    OR.w = OR.w_1 = OR.w_0  = w.thetas
    
    OR.all = list(mu   = OR.mu.part1   + OR.mu.part2,
                  mu_1 = OR.mu_1.part1 + OR.mu_1.part2,
                  mu_0 = OR.mu_0.part1 + OR.mu_0.part2,
                  w   = OR.w,
                  w_1 = OR.w_1,
                  w_0 = OR.w_0)
    
    return(OR.all)
    
  }else{
    
    ### if r > 0: subsampling
    
    # r random binary vectors from A(N.i) with prob dist'n H(.)
    # (N.i x r) matrix, q-th column = a.i.q
    
    a.i.rep <- H.sample(r, X.i, N.i, H.fit)
    
    ### Compute {w/H(a.i.q) + phi/H(a.i.q)}^T F(tau | a.i.q) for q = 1, ... , r
    ### Over mu(Q), mu_1(Q), mu_0(Q) and over thetas and over taus
    ### Resulting matrix `res.taus.thetas` is a list with length = r
    ### res.taus.thetas[q]: list of `mu`, `mu_1`, `mu_0`
    ### res.taus.thetas[q]$`mu`: length(taus) x length(thetas) matrix of 
    ### [{w/H(a.i.q; thetas[d]) + phi/H(a.i.q; thetas[d])}^T F(taus[t] | a.i.q)]_(t,d)
    
    res.thetas.taus <- foreach(q = 1:r) %dopar% OR.over.H.i(a.i.rep[,q])
    
    ### Compute 1/r * sum_{q} {w(a.i.q) + phi(a.i.q)}^T F(tau | a.i.q) / f(a.i.q)
    ### where f(a.i.q) = H(a.i.q,Xi,Ni) 
    ### The result is length(taus) x length(thetas) matrix of 
    ### [1/r * sum_{q} {w/H(a.i.q; thetas[d]) + phi/H(a.i.q; thetas[d])}^T F(taus[t] | a.i.q)]_(t,d)
    
    OR.mu.part2   = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu))
    OR.mu_1.part2 = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu_1))
    OR.mu_0.part2 = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$mu_0))
    
    ### Mean of weights used in subsampling approx, given by
    ### 1/r sum_{q} {w(a.i.q) + phi(a.i.q)}^T 1 / H(a.i.q)
    ### Note: expectation = 1
    OR.w    = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$w  ))
    OR.w_1  = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$w_1))
    OR.w_0  = 1/r * Reduce("+", lapply(res.thetas.taus, function(res) res$w_0))
    
    OR.sub = list(mu   = OR.mu.part1   + OR.mu.part2,
                  mu_1 = OR.mu_1.part1 + OR.mu_1.part2,
                  mu_0 = OR.mu_0.part1 + OR.mu_0.part2,
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

IF.TPB = function(Y.i, D.i, A.i, X.i, N.i, taus, thetas, r = 0, F.fit, G.fit, H.fit){
  
  ### First, compute P(A_bar >= rho | X, N)
  if(r == 0){
    
    ### if r == 0: no subsampling
    # (N.i x 2^N.i) matrix, q-th column = a.i.q
    # P(A_bar >= rho) = \sum_{a.i: a.i.bar >= rho} H(a.i)
    a.i.rep = t(as.matrix(expand.grid(replicate(N.i, 0:1, simplify = FALSE))))
    a.i.bars = colMeans(a.i.rep)
    H.i.rep = apply(X = a.i.rep, MARGIN = 2, FUN = H, X.i = X.i, N.i = N.i, H.fit = H.fit)
    P.A_bar.thetas = sapply(X = thetas, FUN = function(theta) sum(H.i.rep[which(a.i.bars >= theta)]))
    names(P.A_bar.thetas) = thetas
    
  } else {
    
    ### if r > 0: subsampling
    # rr = max(r, 10000) random binary vectors from A(N.i) with prob dist'n H(.)
    # a.i.rep: (N.i x r) matrix, q-th column = a.i.q
    # P(A_bar >= rho) = E[I(A_bar >= rho)] =\approx= 1/rr \sum_{q = 1}^{rr} I(a.i.q_bar >= rho), 
    # where a.i.q \sim H(.)
    
    a.i.rep <- H.sample(max(r,10000), X.i, N.i, H.fit)
    a.i.bars = colMeans(a.i.rep)
    P.A_bar.thetas = sapply(X = thetas, FUN = function(theta) mean(a.i.bars >= theta))
    names(P.A_bar.thetas) = thetas
    
  }
  
  
  ### Wrap up functions for convenience
  F.i.taus = F(taus, A.i, X.i, N.i, F.fit)
  
  ### G.i.Y.i = {G_ij (Y_ij)}_{j=1}^{N_i}
  G.i.Y.i = G.i_Y.i(Y.i, D.i, A.i, X.i, N.i, G.fit)
  
  ### Instead of computing H.i and w.i separately, get w.over.H
  # H.i = H(A.i, X.i, N.i, pi.fit = pi.fit)
  # w.i.thetas = w(A.i, X.i, N.i, thetas, pi.fit = pi.fit)
  w.over.H.i.thetas = w.over.H.TPB(A.i, X.i, N.i, thetas, H.fit, P.A_bar.thetas)
  
  ### martingale integral in AUG term
  mart.i.taus = martingale(Y.i, D.i, A.i, X.i, N.i, taus, F.fit, G.fit)
  
  ### 1. Outcome Regression (OR)
  
  ## IF r >= 2**N.i, then do OR_all instead of OR_sub
  rr = ifelse(r >= 2**N.i, 0, r)
  OR.thetas.taus = OutReg.TPB(A.i, X.i, N.i, taus, thetas, rr, F.fit, H.fit, P.A_bar.thetas)
  
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
