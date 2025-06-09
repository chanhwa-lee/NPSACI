library(dplyr)
library(survival)
library(lme4)
library(randomForestSRC)
library(dbarts)

###---------------------------------###
###---    Training functions     ---###
###---------------------------------###

#' Train a survival model for T (event time) function
#'
#' @param data_train Training dataset with variables Y, D, A, X, N
#' @param X.T.names Character vector of predictor covariate names
#' @param method Model type: either "nonpara" (Random Survival Forest) or "para" (Cox PH)
#' 
#' @return Trained T model object (RSF or Cox PH)
F.train <- function(data_train, X.T.names, method = "nonpara") {
  
  # Construct training data with g.A adjustment
  newdata_train <- data_train %>% 
    group_by(id) %>%
    mutate(g.A = ifelse(n() == 1, 1, (sum(A) - A) / (n() - 1))) %>%
    ungroup() %>%
    select(all_of(c("Y", "D", "A", "g.A", X.T.names)))
  
  # Fit survival model
  if (method == "nonpara") {
    # Nonparametric: Random Survival Forest
    F.fit <- rfsrc(Surv(Y, D) ~ ., 
                   data = newdata_train,
                   perf.type = "none",
                   save.memory = FALSE)
    
  } else if (method == "para") {
    # Parametric: Cox Proportional Hazards
    F.fit <- coxph(Surv(Y, D) ~ ., 
                   data = newdata_train)
    
  } else {
    stop("F.train: Unsupported training method")
  }
  
  return(F.fit)
}


#' Train a survival model for C (censoring) function
#'
#' @param data_train Training dataset with variables Y, D, A, X, N
#' @param X.C.names Character vector of predictor covariate names
#' @param method Model type: either "nonpara" (Random Survival Forest) or "para" (Cox PH)
#' 
#' @return Trained C model object (RSF or Cox PH)
G.train <- function(data_train, X.C.names, method = "nonpara") {
  
  # Construct training data with D <- 1-D modification (RSF cannot take Surv(Y,1-D)) and g.A adjustment
  newdata_train <- data_train %>% 
    mutate(D = 1-D) %>%
    group_by(id) %>%
    mutate(g.A = ifelse(n() == 1, 1, (sum(A) - A) / (n() - 1))) %>%
    ungroup() %>%
    select(all_of(c("Y", "D", "A", "g.A", X.C.names)))
  
  # Fit survival model
  if (method == "nonpara") {
    # Nonparametric: Random Survival Forest
    G.fit <- rfsrc(Surv(Y, D) ~ ., 
                   data = newdata_train,
                   perf.type = "none",
                   save.memory = FALSE)
    
  } else if (method == "para") {
    # Parametric: Cox Proportional Hazards
    G.fit <- coxph(Surv(Y, D) ~ ., 
                   data = newdata_train)
    
  } else {
    stop("G.train: Unsupported training method")
  }
  
  return(G.fit)
}


#' Train a probability model for A (treatment) function
#'
#' @param data_train Training dataset with variables X, N
#' @param X.A.names Character vector of predictor covariate names
#' @param method Model type: either "nonpara" (rBART) or "para" (GLMM)
#' 
#' @return Trained A model object (rBART or GLMM)
H.train <- function(data_train, X.A.names, method = "nonpara") {
  
  # Construct training data
  newdata_train <- data_train %>% 
    select(all_of(c("id", "A", X.A.names)))
  
  # Fit probability model
  if (method == "nonpara") {
    # Nonparametric: rBART
    H.fit = rbart_vi(A ~ . -id, newdata_train, group.by = id,
                     n.samples = 40L, n.burn = 10L, n.thin = 1L,
                     n.chains = 1L,
                     n.trees = 25L, n.threads = 1L)
    
  } else if (method == "para") {
    # Parametric: GLMM
    H.fit = lme4::glmer(A ~ . - id + (1|id), data = newdata_train, family = binomial)
    
  } else {
    stop("H.train: Unsupported training method")
  }
  
  return(H.fit)
}



###---------------------------------###
###---   Prediction functions    ---###
###---------------------------------###


#' Predict the event time function F(t|A,X,N) = P(T ≤ t|A,X,N)
#'
#' @param taus Vector of time points at which to evaluate the survival function
#' @param a.i Treatment vector (length N.i)
#' @param X.i Covariate data frame (N.i x p)
#' @param N.i Number of observations in cluster i
#' @param F.fit Trained T model object (RSF or Cox PH)
#'
#' @return A matrix (N.i x length(taus)) with values P(T ≤ t)
F <- function(taus, a.i, X.i, N.i, F.fit) {
  
  # Create prediction dataset with g.A
  pred_data <- cbind(
    data.frame(A = a.i) %>%
      mutate(g.A = (N.i == 1) * 1 + (N.i != 1) * (sum(A) - A) / (N.i - 1)),
    X.i
  )
  
  # --- Nonparametric estimator: Random Survival Forest ---
  if ("rfsrc" %in% class(F.fit)) {
    S.pred <- predict(F.fit, newdata = pred_data)
    S.times <- S.pred$time.interest
    S.i.S.times <- S.pred$survival
    
    # --- Parametric estimator: Cox PH ---
  } else if ("coxph" %in% class(F.fit)) {
    S.pred <- survfit(F.fit, newdata = pred_data)
    S.times <- S.pred$time
    S.i.S.times <- t(as.matrix(S.pred$surv))
    
  } else {
    stop("F.fit was trained using an unsupported method")
  }
  
  ## S.times: Training times for S.fit
  ## S.i.S.times: (N.i x n.S.time) matrix where (j,k)-element is S_ij(S.times[k]|A_i, X_i, N_i)
  
  # Identify indices of survival times corresponding to taus
  # i.e., for each taus[t], find k s.t. times[k] <= taus[t] < times[k+1]
  indices <- findInterval(taus, S.times)
  
  # Build survival matrix S_i(t): [P(T_ij > taus[t])]_(j,t): (N.i x n.taus) matrix
  S.i.taus = cbind(
    
    ### Potential 0's in the first few elements => P(T_ij > 0) = 1
    matrix(1, nrow = N.i, ncol = sum(indices == 0)),
    
    ### Otherwise, P(T_ij > t) from surv values
    S.i.S.times[,indices[indices != 0]]
  )
  
  ### Return F.i.taus = [P(T_ij <= taus[t])]_(j,t): (N.i x n.taus) matrix
  F.i.taus = 1 - S.i.taus
  colnames(F.i.taus) <- taus
  
  return(F.i.taus)
}


#' Predict the censoring function G(t|A,X,N) = P(C > t | A, X, N)
#'
#' @param taus Vector of time points at which to evaluate the censoring function
#' @param a.i Treatment vector (length N.i)
#' @param X.i Covariate data frame (N.i x p)
#' @param N.i Number of observations in cluster i
#' @param G.fit Trained C model object (RSF or Cox PH)
#'
#' @return A matrix (N.i x length(taus)) with values P(C > t)
G <- function(taus, a.i, X.i, N.i, G.fit) {
  
  # Create prediction dataset with g.A
  pred_data <- cbind(
    data.frame(A = a.i) %>%
      mutate(g.A = (N.i == 1) * 1 + (N.i != 1) * (sum(A) - A) / (N.i - 1)),
    X.i
  )
  
  # --- Nonparametric estimator: Random Survival Forest ---
  if ("rfsrc" %in% class(G.fit)) {
    G.pred <- predict(G.fit, newdata = pred_data)
    G.times <- G.pred$time.interest
    G.i.G.times <- G.pred$survival
    
    # --- Parametric estimator: Cox PH ---
  } else if ("coxph" %in% class(G.fit)) {
    G.pred <- survfit(G.fit, newdata = pred_data)
    G.times <- G.pred$time
    G.i.G.times <- t(as.matrix(G.pred$surv))
    
  } else {
    stop("G.fit was trained using an unsupported method")
  }
  
  ## G.times: Training times for G.fit
  ## G.i.G.times: (N.i x n_G.train) matrix where (j,k)-element is G_ij(G.times[k]|A_i, X_i, N_i)
  
  # Identify indices of survival times corresponding to taus
  # i.e., for each taus[t], find k s.t. times[k] <= taus[t] < times[k+1]
  indices <- findInterval(taus, G.times)
  
  # Build survival matrix G_i(t): [P(C_ij > taus[t])]_(j,t): (N.i x n.taus) matrix
  G.i.taus <- cbind(
    
    ### Potential 0's in the first few elements => P(C_ij > 0) = 1
    matrix(1, nrow = N.i, ncol = sum(indices == 0)),
    
    ### Otherwise, P(C_ij > t) from surv values
    G.i.G.times[, indices[indices != 0]]
  )
  
  ### Return G.i.taus = [P(C_ij > taus[t])]_(j,t): (N.i x n.taus) matrix
  colnames(G.i.taus) <- taus
  
  return(G.i.taus)
}


#' Predict the cluster-level probability H(A|X,N)
#'
#' @param a.i Treatment vector (length N.i)
#' @param X.i Covariate data frame (N.i x p)
#' @param N.i Number of observations in cluster i
#' @param H.fit Trained A model object (RBART or GLMM)
#'
#' @return A scalar representing P(A_i | X_i, N_i)
H <- function(a.i, X.i, N.i, H.fit) {
  
  # --- Nonparametric estimator: Random BART ---
  if ("rbart" %in% class(H.fit)) {
    sigma.b <- mean(H.fit$tau)
    
    f <- function(b) {
      fx.pred <- colMeans(
        predict(
          H.fit, 
          newdata = X.i, 
          group.by = rep(H.fit$group.by[1], N.i), 
          type = "bart"
        )
      )
      pr.b <- pnorm(drop(outer(fx.pred, b, "+")))
      hh <- pr.b^a.i * (1 - pr.b)^(1 - a.i)
      out <- apply(hh, 2, prod) * stats::dnorm(b, mean = 0, sd = sigma.b)
      return(out)
    }
    
    # --- Parametric estimator: Generalized Linear Mixed Model ---
  } else if ("glmerMod" %in% class(H.fit)) {
    beta <- lme4::getME(H.fit, "fixef")
    sigma.b <- lme4::getME(H.fit, "theta")
    
    X <- as.matrix(X.i %>% select(any_of(names(beta))))
    
    f <- function(b) {
      pr.b <- stats::plogis(drop(outer(beta[1] + X %*% beta[-1], b, "+")))
      hh <- pr.b^a.i * (1 - pr.b)^(1 - a.i)
      out <- apply(hh, 2, prod) * stats::dnorm(b, mean = 0, sd = sigma.b)
      return(out)
    }
    
  } else {
    stop("H.fit was trained using an unsupported method")
  }
  
  # Compute integral over b
  H.hat <- stats::integrate(f, lower = -Inf, upper = Inf)$value
  
  return(H.hat)
}



###---------------------------------###
###---     Utility functions     ---###
###---------------------------------###

#' Compute S_i(Y_i) = {S_ij(Y_ij)}_{j=1}^{N_i}
#'
#' @param Y.i Vector of observed times in cluster i (length N.i)
#' @param D.i Vector of event indicators in cluster i (length N.i)
#' @param A.i Vector of treatment assignments in cluster i (length N.i)
#' @param X.i Covariate data frame in cluster i (N.i x p)
#' @param N.i Number of individuals in cluster i
#' @param F.fit Trained T model object (RSF or Cox PH)
#'
#' @return A vector of event time survival probabilities evaluated at Y.i: S_i(Y_i) = {S_ij(Y_ij)}_{j=1}^{N_i}
S.i_Y.i <- function(Y.i, D.i, A.i, X.i, N.i, F.fit) {
  
  # Create prediction dataset with g.A
  pred_data_F <- cbind(
    data.frame(A = A.i) %>%
      mutate(g.A = (N.i == 1) * 1 + (N.i != 1) * (sum(A) - A) / (N.i - 1)),
    X.i
  )
  
  # --- Nonparametric estimator: Random Survival Forest ---
  if ("rfsrc" %in% class(F.fit)) {
    S.pred <- predict(F.fit, newdata = pred_data_F)
    S.times <- S.pred$time.interest
    S.i.S.times <- S.pred$survival
    
    # --- Parametric estimator: Cox PH ---
  } else if ("coxph" %in% class(F.fit)) {
    S.pred <- survfit(F.fit, newdata = pred_data_F)
    S.times <- S.pred$time
    S.i.S.times <- t(as.matrix(S.pred$surv))
    
  } else {
    stop("F.fit was trained using an unsupported method")
  }
  
  # Remove 0 entry in S.i.S.times
  S.i.S.times.col_include_0 = which(apply(S.i.S.times, 2, function(col) any(col == 0)))
  if(length(S.i.S.times.col_include_0) > 0){
    S.times = S.times[-S.i.S.times.col_include_0]
    S.i.S.times = S.i.S.times[,-S.i.S.times.col_include_0]
  }
  
  # --- Initialize output vector ---
  S.i.Y.i <- rep(1, N.i)
  
  # --- Find indices of survival times just before or equal to each Y_ij ---
  indices <- findInterval(Y.i, S.times)                  # N.i-vector
  
  # --- Extract corresponding survival probabilities ---
  for (j in 1:N.i) {
    if (indices[j] != 0) {
      ### times[indices[j]] <= Y_ij < times[indices[j]+1] 
      ### => S_ij(Y_ij) = S_ij(times[indices[j]]) = S.i.S.times[indices[j], j]
      S.i.Y.i[j] <- S.i.S.times[j, indices[j]]
    }
  }
  
  return(S.i.Y.i)
}


#' Compute G_i(Y_i) = {G_ij(Y_ij)}_{j=1}^{N_i}
#'
#' @param Y.i Vector of observed times in cluster i (length N.i)
#' @param D.i Vector of event indicators in cluster i (length N.i)
#' @param A.i Vector of treatment assignments in cluster i (length N.i)
#' @param X.i Covariate data frame in cluster i (N.i x p)
#' @param N.i Number of individuals in cluster i
#' @param G.fit Trained C model object (RSF or Cox PH)
#'
#' @return A vector of censoring survival probabilities evaluated at Y.i: G_i(Y_i) = {G_ij(Y_ij)}_{j=1}^{N_i}
G.i_Y.i <- function(Y.i, D.i, A.i, X.i, N.i, G.fit) {
  
  # Create prediction dataset with g.A
  pred_data_G <- cbind(
    data.frame(A = A.i) %>%
      mutate(g.A = (N.i == 1) * 1 + (N.i != 1) * (sum(A) - A) / (N.i - 1)),
    X.i
  )
  
  # --- Nonparametric estimator: Random Survival Forest ---
  if ("rfsrc" %in% class(G.fit)) {
    G.pred <- predict(G.fit, newdata = pred_data_G)
    G.times <- G.pred$time.interest
    G.i.G.times <- G.pred$survival
    
    # --- Parametric estimator: Cox PH ---
  } else if ("coxph" %in% class(G.fit)) {
    G.pred <- survfit(G.fit, newdata = pred_data_G)
    G.times <- G.pred$time
    G.i.G.times <- t(as.matrix(G.pred$surv))
    
  } else {
    stop("G.fit was trained using an unsupported method")
  }
  
  # Remove 0 entry in G.i.G.times
  G.i.G.times.col_include_0 = which(apply(G.i.G.times, 2, function(col) any(col == 0)))
  if(length(G.i.G.times.col_include_0) > 0){
    G.times = G.times[-G.i.G.times.col_include_0]
    G.i.G.times = G.i.G.times[,-G.i.G.times.col_include_0]
  }
  
  # --- Initialize output vector ---
  G.i.Y.i <- rep(1, N.i)
  
  # --- Find indices of survival times just before or equal to each Y_ij ---
  indices <- findInterval(Y.i, G.times)
  
  # --- Extract corresponding survival probabilities ---
  for (j in 1:N.i) {
    if (indices[j] != 0) {
      ### times[indices[j]] <= Y_ij < times[indices[j]+1] 
      ### => G_ij(Y_ij) = G_ij(times[indices[j]]) = G.i.G.times[indices[j], j]
      G.i.Y.i[j] <- G.i.G.times[j, indices[j]]
    }
  }
  
  return(G.i.Y.i)
}


#' Compute Martingale integral in AUG term
#'
#' @param Y.i Vector of observed times in cluster i (length N.i)
#' @param D.i Vector of event indicators in cluster i (length N.i)
#' @param A.i Vector of treatment assignments in cluster i (length N.i)
#' @param X.i Covariate data frame in cluster i (N.i x p)
#' @param N.i Number of individuals in cluster i
#' @param taus Vector of time points at which to evaluate the survival function 
#' @param F.fit Trained T model object (RSF or Cox PH)
#' @param G.fit Trained C model object (RSF or Cox PH)
#'
#' @return A vector of stieltjes integrals in AUG terms

martingale = function(Y.i, D.i, A.i, X.i, N.i, taus, F.fit, G.fit){
  
  #--- T model ingredients ---
  pred_data_F <- cbind(
    data.frame(A = A.i) %>%
      mutate(g.A = (N.i == 1) * 1 + (N.i != 1) * (sum(A) - A) / (N.i - 1)),
    X.i
  )
  
  # --- Nonparametric estimator: Random Survival Forest ---
  if ("rfsrc" %in% class(F.fit)) {
    S.pred <- predict(F.fit, newdata = pred_data_F)
    S.times <- S.pred$time.interest
    S.i.S.times <- S.pred$survival
    
    # --- Parametric estimator: Cox PH ---
  } else if ("coxph" %in% class(F.fit)) {
    S.pred <- survfit(F.fit, newdata = pred_data_F)
    S.times <- S.pred$time
    S.i.S.times <- t(as.matrix(S.pred$surv))
    
  } else {
    stop("F.fit was trained using an unsupported method")
  }
  
  # Remove 0 entry in S.i.S.times
  S.i.S.times.col_include_0 = which(apply(S.i.S.times, 2, function(col) any(col == 0)))
  if(length(S.i.S.times.col_include_0) > 0){
    S.times = S.times[-S.i.S.times.col_include_0]
    S.i.S.times = S.i.S.times[,-S.i.S.times.col_include_0]
  }
  
  
  #--- C model ingredients ---
  pred_data_G <- cbind(
    data.frame(A = A.i) %>%
      mutate(g.A = (N.i == 1) * 1 + (N.i != 1) * (sum(A) - A) / (N.i - 1)),
    X.i
  )
  
  
  # --- Nonparametric estimator: Random Survival Forest ---
  if ("rfsrc" %in% class(G.fit)) {
    G.pred <- predict(G.fit, newdata = pred_data_G)
    G.times <- G.pred$time.interest
    G.i.G.times <- G.pred$survival
    
    # --- Parametric estimator: Cox PH ---
  } else if ("coxph" %in% class(G.fit)) {
    G.pred <- survfit(G.fit, newdata = pred_data_G)
    G.times <- G.pred$time
    G.i.G.times <- t(as.matrix(G.pred$surv))
    
  } else {
    stop("G.fit was trained using an unsupported method")
  }
  
  # Remove 0 entry in G.i.G.times
  G.i.G.times.col_include_0 = which(apply(G.i.G.times, 2, function(col) any(col == 0)))
  if(length(G.i.G.times.col_include_0) > 0){
    G.times = G.times[-G.i.G.times.col_include_0]
    G.i.G.times = G.i.G.times[,-G.i.G.times.col_include_0]
  }
  
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


#' Sample binary vectors from the cluster-level probability H(.|X,N)
#'
#' @param r Number of binary vectors to be sampled
#' @param X.i Covariate data frame (N.i x p)
#' @param N.i Number of observations in cluster i
#' @param H.fit Trained A model object (RBART or GLMM)
#'
#' @return A matrix of sampled binary vectors (N.i x r)
H.sample = function(r, X.i, N.i, H.fit){
  
  # r random binary vectors from A(N.i) with prob dist'n H(.|X.i,N.i)
  # a.i.q \sim H(.) <=> b.i.q \sim N(0, sigma^2), a.ij.q|b.i.q \sim Ber(pi(X.ij, b.i.q))
  # Return (N.i x r) matrix, q-th column = a.i.q
  
  # --- Nonparametric estimator: Random BART ---
  if ("rbart" %in% class(H.fit)) {
    sigma.b <- mean(H.fit$tau)
    b.rep = rnorm(r, 0, sigma.b)
    fx.pred = colMeans(predict(H.fit, newdata = X.i, group.by = rep(H.fit$group.by[1], N.i), type = "bart"))
    pr.b.rep <- pnorm(drop(outer(fx.pred, b.rep, "+")))
    
    # --- Parametric estimator: Generalized Linear Mixed Model ---
  } else if ("glmerMod" %in% class(H.fit)) {
    sigma.b = lme4::getME(H.fit, "theta")
    b.rep = rnorm(r, 0, sigma.b)
    beta = lme4::getME(H.fit, "fixef")
    X = as.matrix(X.i %>% select(any_of(names(beta))))
    fx.pred = beta[1] + X %*% beta[-1]
    pr.b.rep <- stats::plogis(drop(outer(fx.pred, b.rep, "+")))
    
  } else {
    stop("H.fit was trained using an unsupported method")
  }
  
  a.i.rep <- apply(pr.b.rep, 2, function(p) rbinom(N.i, 1, p))
  return(a.i.rep)
  
}
