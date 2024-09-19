library(survival)
library(dplyr)
library(lme4)

###---------------------------------###
###---     Survival functions    ---###
###---------------------------------###

###--- F (survival) function training ---###
#' @input Training data O_train = {(Y, D, A, X, N)}_train
#' @intermediate 
#' @output survival function fit object
#' @example F.fit = F.train(data_train, X.T.names)
F.train = function(data_train, X.T.names){
  
  ## Fit Cox PH regression
  F.fit <- coxph(Surv(Y,D) ~ ., 
                 data = data_train %>% group_by(id) %>% 
                   mutate(g.A = (n()==1)*1 + (n()!=1)*(sum(A) - A)/(n()-1)) %>% ungroup() %>%
                   select(c("Y", "D", "A", "g.A", X.T.names)))
  
  return(F.fit)
}


###--- F (survival) function prediction ---###
#' @input Triplet (a.i, X.i, N.i), times of interest taus
#' @intermediate F.train object
#' @output survival function prediction
#' @example F(taus, a.i, X.i, N.i)
F = function(taus, a.i, X.i, N.i, F.fit){
  
  S.pred = survfit(F.fit, 
                   newdata = cbind(data.frame(A = a.i) %>% 
                                     mutate(g.A = (N.i==1)*1 + (N.i!=1)*(sum(A) - A)/(N.i-1)), X.i))
  
  ### Training times for S.fit
  S.times = S.pred$time
  
  ## S.i.S.times: (N.i x n.S.time) matrix where (j,k)-element is S_ij(S.times[k]|A_i, X_i, N_i)
  S.i.S.times = t(as.matrix(S.pred$surv))
  
  ### Find the first element in times that is less than or equal to each element in taus
  ### i.e., for each taus[t], find k s.t. times[k] <= taus[t] < times[k+1]
  indices <- findInterval(taus, S.times)
  
  ### Fitted [P(T_ij > taus[t])]_(j,t): (N.i x n.taus) matrix
  S.i.taus = cbind(
    
    ### Potential 0's in the first few elements => P(T_ij > 0) = 1
    matrix(1, nrow = N.i, ncol = sum(indices == 0)),
    
    ### Otherwise, P(T_ij > t) from surv values
    S.i.S.times[,indices[indices != 0]]
  )
  
  ### Return F.i.taus = [P(T_ij <= taus[t])]_(j,t): (N.i x n.taus) matrix
  F.i.taus = 1 - S.i.taus
  
  colnames(F.i.taus) = taus
  
  return(F.i.taus)
  
}


###--- G (censoring) function fitting ---###
#' @input Training data O_train = {(Y, D, A, X, N)}_train
#' @intermediate 
#' @output censoring function fit object and training time points (for stieltjes integral)
#' @example 
#' G.fit = G.train(data_train, X.C.names)
G.train = function(data_train, X.C.names){
  
  ## Fit Cox PH regression
  G.fit <- coxph(Surv(Y,1-D) ~ ., 
                 data = data_train %>% select(c("Y","D","A", X.C.names)))
  
  return(G.fit)
}


###--- G (censoring) function prediction ---###
#' @input Triplet (a.i, X.i, N.i), times of interest taus
#' @intermediate G.train object
#' @output censoring function prediction 
G = function(taus, a.i, X.i, N.i, G.fit){
  
  G.pred = survfit(G.fit, 
                   newdata = cbind(data.frame(A = a.i), X.i))
  
  ### Training times for G.fit
  G.times = G.pred$time
  
  ## G.i.G.times: (N.i x n_G.train) matrix where (j,k)-element is G_ij(G.times[k]|A_i, X_i, N_i)
  G.i.G.times = t(as.matrix(G.pred$surv))
  
  ### Find the first element in times that is less than or equal to each element in taus
  ### i.e., for each taus[t], find k s.t. times[k] <= taus[t] < times[k+1]
  indices <- findInterval(taus, G.times)
  
  ### Fitted [P(C_ij > taus[t])]_(j,t): (N.i x n.taus) matrix
  G.i.taus = cbind(
    
    ### Potential 0's in the first few elements => P(C_ij > 0) = 1
    matrix(1, nrow = N.i, ncol = sum(indices == 0)),
    
    ### Otherwise, P(C_ij > t) from surv values
    G.i.G.times[, indices[indices != 0]]
  )
  
  ### Return G.i.G.times_at_taus = [P(C_ij > taus[t])]_(j,t): (N.i x n.taus) matrix
  colnames(G.i.taus) = taus
  
  return(G.i.taus)
  
}


###--- S.i(Y.i) function ---###
#' @input cluster level data (Y.i, D.i, A.i, X.i, N.i), 
#' nuisance estimators F.fit, 
#' @intermediate F, training event times
#' @output S.i(Y.i) = {S_ij(Y_ij)}_{j=1}^{N_i}

S.i_Y.i = function(Y.i, D.i, A.i, X.i, N.i, F.fit){
  
  ### Prediction over new data from F.fit
  S.pred = survfit(F.fit, 
                   newdata = cbind(data.frame(A = A.i) %>% 
                                     mutate(g.A = (N.i==1)*1 + (N.i!=1)*(sum(A) - A)/(N.i-1)), X.i))
  
  ### Training times for S.fit
  S.times = S.pred$time
  
  ## S.i.S.times: (N.i x n.S.time) matrix where (j,k)-element is S_ij(S.times[k]|A_i, X_i, N_i)
  S.i.S.times = t(as.matrix(S.pred$surv))
  
  ### Find the first element in times that is less than or equal to each element in Y.i
  ### i.e., for each Y_ij, find k s.t. S.times[k] <= Y_ij = Y.i[j] < S.times[k+1]
  indices <- findInterval(Y.i, S.times) ## N.i-vector
  
  ### S.i.Y.i = {S_ij (Y_ij)}_{j=1}^{N_i} N.i-vector
  S.i.Y.i <- rep(1, N.i)
  for(j in 1:N.i){
    if(indices[j] != 0){
      ### times[indices[j]] <= Y_ij < times[indices[j]+1] 
      ### => S_ij(Y_ij) = S_ij(times[indices[j]]) = S.i.S.times[indices[j], j]
      S.i.Y.i[j] = S.i.S.times[j, indices[j]]
    }
  }
  
  return(S.i.Y.i)
  
}

###--- G.i(Y.i) function ---###
#' @input cluster level data (Y.i, D.i, A.i, X.i, N.i), 
#' nuisance estimators G.fit, 
#' @intermediate G, training event times
#' @output G.i(Y.i) = {G_ij(Y_ij)}_{j=1}^{N_i}

G.i_Y.i = function(Y.i, D.i, A.i, X.i, N.i, G.fit){
  
  ### Prediction over new data from G.fit
  G.pred = survfit(G.fit, 
                   newdata = cbind(data.frame(A = A.i), X.i))
  
  ### Training times for G.fit
  G.times = G.pred$time
  
  ## G.i.G.times: (N.i x n_G.train) matrix where (j,k)-element is G_ij(G.times[k]|A_i, X_i, N_i)
  G.i.G.times = t(as.matrix(G.pred$surv))
  
  ### Find the first element in times that is less than or equal to each element in Y.i
  ### i.e., for each Y_ij, find k s.t. G.times[k] <= Y_ij = Y.i[j] < G.times[k+1]
  indices <- findInterval(Y.i, G.times) ## N.i-vector
  
  ### G.i.Y.i = {G_ij (Y_ij)}_{j=1}^{N_i}: N.i-vector
  G.i.Y.i <- rep(1, N.i)
  for(j in 1:N.i){
    if(indices[j] != 0){
      ### times[indices[j]] <= Y_ij < times[indices[j]+1] 
      ### => G_ij(Y_ij) = G_ij(times[indices[j]]) = G.i.G.times[indices[j], j]
      G.i.Y.i[j] = G.i.G.times[j,indices[j]]
    }
  }
  
  return(G.i.Y.i)
  
}


###--- Martingale integral function ---###
#' @input cluster level data (Y.i, D.i, A.i, X.i, N.i), 
#' nuisance estimators (F.fit, G.fit), 
#' times of interest \taus
#' @intermediate F, G, training event times
#' @output stieltjes integral in augmented term for individual j in cluster i

martingale = function(Y.i, D.i, A.i, X.i, N.i, taus, F.fit, G.fit){
  
  ### Prediction over new data from F.fit
  S.pred = survfit(F.fit, 
                   newdata = cbind(data.frame(A = A.i) %>% 
                                     mutate(g.A = (N.i==1)*1 + (N.i!=1)*(sum(A) - A)/(N.i-1)), X.i))
  
  ### Training times for S.fit
  S.times = S.pred$time
  
  ## S.i.S.times: (N.i x n.S.time) matrix where (j,k)-element is S_ij(S.times[k]|A_i, X_i, N_i)
  S.i.S.times = t(as.matrix(S.pred$surv))
  

  ### Prediction over new data from G.fit
  G.pred = survfit(G.fit, 
                   newdata = cbind(data.frame(A = A.i), X.i))
  
  ### Training times for G.fit
  G.times = G.pred$time
  
  ## G.i.G.times: (N.i x n_G.train) matrix where (j,k)-element is G_ij(G.times[k]|A_i, X_i, N_i)
  G.i.G.times = t(as.matrix(G.pred$surv))
  
  #### First part: Integral wrt censoring process (N_ij^C(r))
  
  ### part1: (N.i x n.taus) matrix
  ### [part1]_(j,t) = (1-Delta_ij) * 1/G_ij(Y_ij) * (1-S_ij(taus[t])/S_ij(Y_ij)) *I(Y_ij <= taus[t])
  ### j = 1, ..., N.i; t = 1, ..., n.taus
  
  ### S.i.Y.i = {S_ij (Y_ij)}_{j=1}^{N_i}
  ### G.i.Y.i = {G_ij (Y_ij)}_{j=1}^{N_i}
  S.i.Y.i = S.i_Y.i(Y.i, D.i, A.i, X.i, N.i, F.fit)
  G.i.Y.i = G.i_Y.i(Y.i, D.i, A.i, X.i, N.i, G.fit)
  
  ### (N.i x n.taus) matrix where (j,t)-element is S_ij(taus[t])
  S.i.taus = 1 - F(taus, A.i, X.i, N.i, F.fit)
  
  ### (N.i x n.taus) matrix
  part1 =
    (1 - D.i) *
    1/G.i.Y.i *
    (1-S.i.taus / S.i.Y.i) *
    outer(Y.i, taus,  "<=")           ## [I(Y_ij <= taus[t])]_(j,t) (N.i x n.taus) matrix
  
  colnames(part1) = taus
  
  ### CHECK: part1 is always non-negative and increasing in taus
  
  #### Second part: Stieltjes integral, summation at jumps, which are G.fit training times
  #### Idea: decompose part2 into multiplication of matrices
  
  ### part2: (N.i x n.taus) matrix
  ### [part2]_(j,t) = \sum_{k= 1,..., n.G.times}
  ###                      1/G_ij(G.times[k])^2
  ###                      * {G_ij(G.times[k]) - G_ij(G.times[k-1]}
  ###                      * 1(G.times[k] <= Y_ij)
  ###                      * 1(G.times[k] <= taus[t])
  ###                      * (1 - S_ij(taus[t]) / S_ij(G.times[k]))
  ### j = 1, ..., N.i; t = 1, ..., n.taus
  
  ### Decompose: part2 = part2_1 - part2_2
  
  ### [part2_1]_(j,t) = \sum_{k= 1,..., n.G.times}
  ###                      1/G_ij(G.times[k])^2
  ###                      * {G_ij(G.times[k]) - G_ij(G.times[k-1]}
  ###                      * 1(G.times[k] <= Y_ij)
  ###                      * 1(G.times[k] <= taus[t])
  
  ### [part2_2]_(j,t) = \sum_{k= 1,..., n.G.times}
  ###                      1/G_ij(G.times[k])^2
  ###                      * {G_ij(G.times[k]) - G_ij(G.times[k-1]}
  ###                      * 1(G.times[k] <= Y_ij)
  ###                      * 1(G.times[k] <= taus[t])
  ###                      * 1 / S_ij(G.times[k])
  ###                      * S_ij(taus[t])
  
  ### Observation: Each term is summation over k. 
  
  ### part2_1: can be expressed by
  ###   c_jt = \sum_k a_jk b_kt
  ### Then,
  ### C = [c_jt]_(j,t) = [a_jk]_(j,k) %*% [b_kt]_(k,t) = A %*% B (%*%: matrix multiplication)
  
  ### a_jk = 1/G_ij(G.times[k])^2  <---------------------- G.i.G.times
  ###        * {G_ij(G.times[k]) - G_ij(G.times[k-1]}  <-- G.i.G.times.diff
  ###        * 1(G.times[k] <= Y_ij)  <------------------- I.Y.i.ge.G.times
  ### b_kt = 1(G.times[k] <= taus[t]) <------------------- I.G.times.le.taus
  
  G.i.G.times.diff = G.i.G.times - cbind(1, G.i.G.times[,-length(G.times)])
  I.Y.i.ge.G.times = outer(Y.i, G.times, ">=")
  I.G.times.le.taus = outer(G.times, taus, "<=")
  
  part2_1 = (1/G.i.G.times^2 * G.i.G.times.diff * I.Y.i.ge.G.times) %*% I.G.times.le.taus
  
  ### CHECK: part2_1 should be decreasing over tau at each row and always negative. 
  ###        Furthermore, should be constant once Y.ij < taus[t]
  
  ### part2_2: can be expressed by
  ###   e_jt = \sum_k a'_jk b_kt d_jt
  ### Then,
  ### E = [e_jt]_(j,t) = {[a'_jk]_(j,k) %*% [b_kt]_(k,t)} * [d_jt]_(j,t) = (A' %*% B) * D (*: element-wise multiplication)
  
  ### a'_jk = a_jk
  ###        * 1 / S_ij(G.times[k])  <--------------------- S.i.G.times
  ### b_kt = 1(G.times[k] <= taus[t])  <------------------- I.G.times.le.taus
  ### d_jt = S_ij(taus[t])  <------------------------------ S.i.taus
  
  ### Need to compute S.i.G.times = [S_ij(G.times[k])]_(j,k): (N.i x n.G.times) matrix
  indices <- findInterval(G.times, S.times) ## (n.G.times)-vector
  
  S.i.G.times = cbind(
    
    ### Potential 0's in the first few elements => P(T_ij > 0) = 1
    matrix(1, nrow = N.i, ncol = sum(indices == 0)),
    
    ### Otherwise, P(T_ij > t) from surv values
    S.i.S.times[, indices[indices != 0]]
  )
  
  part2_2 = ((1/G.i.G.times^2 * G.i.G.times.diff * I.Y.i.ge.G.times * 1/S.i.G.times) %*% I.G.times.le.taus) * S.i.taus
  
  ### CHECK: part2_2 should be always negative.
  ###        Once Y.ij < taus[t], decreasing in magnitude. Until Y.ij > taus[t], no clear trend
  
  part2 = part2_1 - part2_2
  
  colnames(part2) = taus
  
  ### CHECK: part2 should be always negative and decrease in taus
  
  #### Return sum of two parts
  #### (N.i x n.taus) matrix
  return(part1 + part2)
  
}


###---------------------------------###
###---     Treatment function    ---###
###---------------------------------###

###--- H (cluster propensity score) function fitting ---###
#' @input Training data O_train = {(A, X, N)}_train
#' @output cluster probability function fit object
H.train = function(data_train, X.A.names){
  
  H.fit = lme4::glmer(A ~ . - id + (1|id), data = data_train %>% select(c("id", "A", X.A.names)), family = binomial)
  
  return(H.fit)
}

###--- H (cluster probability) function prediction ---###
#' @input Triplet (a.i, X.i, N.i)
#' @intermediate possibly pi (individual-level propensity score)
#' @output
H = function(a.i, X.i, N.i, H.fit){

  beta = lme4::getME(H.fit, "fixef")
  sigma.b = lme4::getME(H.fit, "theta")
  
  X = as.matrix(X.i %>% select(any_of(names(beta))))
  
  f = function(b){
    pr.b <- stats::plogis(drop(outer(beta[1] + X %*% beta[-1], b, "+")))
    hh <- (pr.b)^a.i * (1 - pr.b)^(1 - a.i)
    out <- apply(hh, 2, prod) * stats::dnorm(b, mean = 0, sd = sigma.b)
    return(out)
  }
  
  H.hat = stats::integrate(f, lower = -Inf, upper = Inf)$value
  
  return(H.hat)
  
}



###--- Sampling binary vectors from H (cluster probability) ---###
#' @input Triplet (a.i, X.i, N.i)
#' @intermediate possibly pi (individual-level propensity score)
#' @output
H.sample = function(r, X.i, N.i, H.fit){
  
  # r random binary vectors from A(N.i) with prob dist'n H(.)
  # a.i.q \sim H(.) <=> b.i.q \sim N(0, sigma^2), a.ij.q|b.i.q \sim Ber(expit(X.ij beta + b.i.q))
  # (N.i x r) matrix, q-th column = a.i.q
    
  sigma.b = lme4::getME(H.fit, "theta")
  
  b.rep = rnorm(r, 0, sigma.b)
  
  beta = lme4::getME(H.fit, "fixef")
  X = as.matrix(X.i %>% select(any_of(names(beta))))
  
  fx.pred = beta[1] + X %*% beta[-1]
  
  pr.b.rep <- stats::plogis(drop(outer(fx.pred, b.rep, "+")))
  
  a.i.rep <- apply(pr.b.rep, 2, function(p) rbinom(N.i, 1, p))
  
  return(a.i.rep)
  
}
