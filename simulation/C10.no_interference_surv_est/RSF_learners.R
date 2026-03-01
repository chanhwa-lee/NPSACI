#------------------------------------------------------------------------------#
# DATA STRUCTURE REQUIREMENTS
# 1. "Time": Numeric (Observed time)
# 2. "Event": Binary (0=Censored, 1=Event)
# 3. "Treatment": Binary (0=Control, 1=Treated)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# 0. Libraries & Setup
#------------------------------------------------------------------------------#
library(randomForestSRC) 
library(pec)             
library(gbm)             
library(randomForest)    
library(riskRegression)  
library(survival)        
library(grf)             

#------------------------------------------------------------------------------#
# 1. T-Learner with Bootstrap CI
#------------------------------------------------------------------------------#
get_ATE_T_learner <- function(data, testdat=NULL, time.interest, n.boot=20){
  # Default testdat
  if(is.null(testdat)) testdat <- data
  data <- data[complete.cases(data),]
  
  # Helper function to calculate ATE once
  calc_ate <- function(curr_data, curr_test){
    data1 <- curr_data[curr_data$Treatment==1, -which(names(curr_data)=="Treatment")]
    data0 <- curr_data[curr_data$Treatment==0, -which(names(curr_data)=="Treatment")]
    
    rfsrc_1 <- rfsrc(Surv(Time, Event) ~ ., data = data1)
    rfsrc_0 <- rfsrc(Surv(Time, Event) ~ ., data = data0)
    
    p1 <- predictSurvProb(rfsrc_1, newdata = curr_test, times=time.interest)
    p0 <- predictSurvProb(rfsrc_0, newdata = curr_test, times=time.interest)
    return(mean(p1 - p0, na.rm=TRUE))
  }
  
  # 1. Point Estimate
  point_est <- calc_ate(data, testdat)
  
  # 2. Bootstrap for CI
  boot_ests <- numeric(n.boot)
  n <- nrow(data)
  
  for(i in 1:n.boot){
    # Resample data with replacement
    idx <- sample(1:n, n, replace=TRUE)
    boot_data <- data[idx, ]
    # We predict on the original testdat to see variation in the model
    boot_ests[i] <- calc_ate(boot_data, testdat)
  }
  
  ci <- quantile(boot_ests, probs = c(0.025, 0.975))
  
  return(list(ate = point_est, ci_lower = ci[1], ci_upper = ci[2]))
}

#------------------------------------------------------------------------------#
# 2. X-Learner with Bootstrap CI
#------------------------------------------------------------------------------#
get_ATE_X_learner <- function(data, testdat=NULL, time.interest, ntrees=500, n.boot=20){
  # Default testdat
  if(is.null(testdat)) testdat <- data
  data <- data[complete.cases(data),]
  
  # Helper function
  calc_ate <- function(curr_data, curr_test){
    data1 <- curr_data[curr_data$Treatment==1, -which(names(curr_data)=="Treatment")]
    data0 <- curr_data[curr_data$Treatment==0, -which(names(curr_data)=="Treatment")]
    
    # Base learners
    rf1 <- rfsrc(Surv(Time, Event) ~ ., data = data1)
    rf0 <- rfsrc(Surv(Time, Event) ~ ., data = data0)
    
    # Imputed outcomes
    p0_on_1 <- predictSurvProb(rf0, newdata = data1, times=time.interest)
    p1_on_1 <- predictSurvProb(rf1, newdata = data1, times=time.interest)
    data1$d1 <- p1_on_1 - p0_on_1
    
    p0_on_0 <- predictSurvProb(rf0, newdata = data0, times=time.interest)
    p1_on_0 <- predictSurvProb(rf1, newdata = data0, times=time.interest)
    data0$d0 <- p1_on_0 - p0_on_0
    
    # GBMs
    cols_ex <- c("Time", "Event", "sim", "status", "Treatment")
    m0 <- gbm(d0 ~ ., data = data0[!is.na(data0$d0), !names(data0) %in% cols_ex], 
              distribution="gaussian", n.trees = ntrees, interaction.depth=4, n.minobsinnode = 5, bag.fraction = 0.8)
    m1 <- gbm(d1 ~ ., data = data1[!is.na(data1$d1), !names(data1) %in% cols_ex], 
              distribution="gaussian", n.trees = ntrees, interaction.depth=4, n.minobsinnode = 5, bag.fraction = 0.8)
    
    pred_d0 <- predict(m0, curr_test, n.trees=ntrees)
    pred_d1 <- predict(m1, curr_test, n.trees=ntrees)
    
    # Propensity
    d_prop <- curr_data[, !names(curr_data) %in% c("d0", "d1")]
    d_prop$Treatment <- as.factor(d_prop$Treatment)
    rf_p <- randomForest(Treatment ~ ., data=d_prop, ntree = 100) # Faster RF for prop
    g_x <- predict(rf_p, curr_test, type="prob")[,2]
    
    cate <- g_x * pred_d0 + (1 - g_x) * pred_d1
    return(mean(cate, na.rm=TRUE))
  }
  
  # 1. Point Estimate
  point_est <- calc_ate(data, testdat)
  
  # 2. Bootstrap
  boot_ests <- numeric(n.boot)
  n <- nrow(data)
  
  for(i in 1:n.boot){
    idx <- sample(1:n, n, replace=TRUE)
    boot_data <- data[idx, ]
    boot_ests[i] <- calc_ate(boot_data, testdat)
  }
  
  ci <- quantile(boot_ests, probs = c(0.025, 0.975))
  
  return(list(ate = point_est, ci_lower = ci[1], ci_upper = ci[2]))
}

#------------------------------------------------------------------------------#
# 3. AIPTW (Analytical CI)
#------------------------------------------------------------------------------#
get_ATE_AIPTW <- function(data, time.interest){
  data$Treatment <- as.factor(data$Treatment)
  
  ate_fit <- ate(
    event = Surv(Time, Event) ~ ., 
    treatment = Treatment ~ ., 
    censor = Surv(Time, Event == 0) ~ .,
    data = data,
    times = time.interest,
    cause = 1,
    estimator = "AIPTW",
    cause.1.nuisance.handler = "rfsrc", 
    propensity.nuisance.handler = "forest",
    censoring.nuisance.handler = "rfsrc",
    verbose = FALSE
  )
  
  # ate_fit calculates RISK (Probability of Event).
  # We want SURVIVAL (Probability of No Event).
  # Survival = 1 - Risk.
  # Difference(Surv) = (1 - R1) - (1 - R0) = R0 - R1 = -(R1 - R0).
  # Confidence Interval flips: [-Upper_Risk, -Lower_Risk]
  
  risk_diff_est <- ate_fit$diffRisk$estimate
  risk_diff_lower <- ate_fit$diffRisk$lower
  risk_diff_upper <- ate_fit$diffRisk$upper
  
  return(list(
    ate = -risk_diff_est, 
    ci_lower = -risk_diff_upper, 
    ci_upper = -risk_diff_lower
  ))
}

#------------------------------------------------------------------------------#
# 4. CSF (Analytical CI using GRF)
#------------------------------------------------------------------------------#
get_ATE_CSF <- function(data, testdat=NULL, time.interest){
  if(is.null(testdat)) testdat <- data
  
  cols_exclude <- c("Time", "Event", "Treatment")
  X_train <- as.matrix(data[, !names(data) %in% cols_exclude])
  Y_train <- data$Time
  D_train <- data$Event
  W_train <- as.numeric(as.character(data$Treatment))
  
  cs_forest <- causal_survival_forest(
    X = X_train, Y = Y_train, D = D_train, W = W_train,
    target = "survival.probability", horizon = time.interest
  )
  
  # Calculate Average Treatment Effect on the TRAINING sample (Population ATE)
  # GRF provides efficient standard errors via 'average_treatment_effect'
  ate_obj <- average_treatment_effect(cs_forest, target.sample = "all")
  
  est <- ate_obj["estimate"]
  se  <- ate_obj["std.err"]
  
  return(list(
    ate = est, 
    ci_lower = est - 1.96 * se, 
    ci_upper = est + 1.96 * se
  ))
}

#------------------------------------------------------------------------------#
# 5. MISTR (Rubin's Rules for CI)
#------------------------------------------------------------------------------#
get_ATE_MISTR_R <- function(data, testdat=NULL, time.interest, A=5, n.trees=500){
  if(is.null(testdat)) testdat <- data
  
  # 1. Imputation Model
  rf_impute <- rfsrc(Surv(Time, Event) ~ ., data = data, ntree = n.trees)
  censored_idx <- which(data$Event == 0 & data$Time < time.interest)
  surv_obj <- predictSurvProb(rf_impute, newdata = data[censored_idx,], times = rf_impute$time.interest)
  model_times <- rf_impute$time.interest
  
  # Vectors to store ATE and Variance for each imputation
  ate_within <- numeric(A)
  var_within <- numeric(A)
  
  for(a in 1:A){
    imputed_data <- data
    
    # Imputation Loop
    for(i in seq_along(censored_idx)){
      idx <- censored_idx[i]
      C_i <- data$Time[idx] 
      valid_time_indices <- which(model_times > C_i)
      
      if(length(valid_time_indices) > 0){
        closest_t_idx <- which.min(abs(model_times - C_i))
        S_C <- surv_obj[i, closest_t_idx]
        if(S_C < 1e-6) S_C <- 1e-6 
        cond_probs <- surv_obj[i, valid_time_indices] / S_C
        
        u <- runif(1)
        drop_idx <- which(cond_probs < u)[1]
        
        if(!is.na(drop_idx)){
          v_time <- model_times[valid_time_indices[drop_idx]]
        } else {
          v_time <- max(model_times) + 0.1
        }
        imputed_data$Time[idx]  <- v_time
        imputed_data$Event[idx] <- 1 
      }
    }
    
    # Estimation (Causal Forest)
    Y_bin <- as.numeric(imputed_data$Time > time.interest)
    cols_exclude <- c("Time", "Event", "Treatment")
    X_mat <- as.matrix(imputed_data[, !names(imputed_data) %in% cols_exclude])
    W_vec <- as.numeric(as.character(imputed_data$Treatment)) 
    
    cf <- causal_forest(X = X_mat, Y = Y_bin, W = W_vec)
    
    # Get ATE and Variance for THIS imputation using GRF's built-in inference
    res <- average_treatment_effect(cf)
    ate_within[a] <- res["estimate"]
    var_within[a] <- res["std.err"]^2
  }
  
  # RUBIN'S RULES
  # 1. Pooled Estimate (Mean of estimates)
  ate_pooled <- mean(ate_within)
  
  # 2. Within Variance (Mean of variances)
  mean_var_within <- mean(var_within)
  
  # 3. Between Variance (Variance of estimates)
  var_between <- var(ate_within)
  
  # 4. Total Variance
  var_total <- mean_var_within + (1 + 1/A) * var_between
  
  # 5. Total Standard Error
  se_total <- sqrt(var_total)
  
  return(list(
    ate = ate_pooled,
    ci_lower = ate_pooled - 1.96 * se_total,
    ci_upper = ate_pooled + 1.96 * se_total
  ))
}


#------------------------------------------------------------------------------#
# Main Wrapper Function
#------------------------------------------------------------------------------#
get_all_ATE_results <- function(data, time.interest, n.boot=20, A_mistr=5){
  
  print("--- Starting Comprehensive ATE Analysis ---")
  
  methods_list <- list()
  
  # Helper to track time and print it
  run_method <- function(name, expr) {
    print(paste0("Running ", name, "..."))
    start_time <- Sys.time()
    result <- expr
    end_time <- Sys.time()
    
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    print(paste0("   -> Finished in ", round(duration, 2), " seconds."))
    
    return(result)
  }
  
  # 1. T-Learner
  methods_list[["T-Learner"]] <- run_method("1. T-Learner (Bootstrap CI)", 
                                            get_ATE_T_learner(data, time.interest=time.interest, n.boot=n.boot))
  
  # 2. X-Learner
  methods_list[["X-Learner"]] <- run_method("2. X-Learner (Bootstrap CI)", 
                                            get_ATE_X_learner(data, time.interest=time.interest, n.boot=n.boot))
  
  # 3. AIPTW
  methods_list[["AIPTW"]] <- run_method("3. AIPTW (analytical CI)", 
                                        get_ATE_AIPTW(data, time.interest=time.interest))
  
  # 4. CSF
  methods_list[["CSF"]] <- run_method("4. CSF (EIF CI)", 
                                      get_ATE_CSF(data, time.interest=time.interest))
  
  # 5. MISTR
  methods_list[["MISTR"]] <- run_method("5. MISTR (Rubin's Rule CI)", 
                                        get_ATE_MISTR_R(data, time.interest=time.interest, A=A_mistr))
  
  # Combine into Dataframe
  results_df <- do.call(rbind, lapply(names(methods_list), function(m) {
    r <- methods_list[[m]]
    data.frame(
      Method = m,
      Estimate = r$ate,
      Lower_95_CI = r$ci_lower,
      Upper_95_CI = r$ci_upper,
      stringsAsFactors = FALSE
    )
  }))
  
  # Rounding for cleanliness
  results_df$Estimate    <- round(results_df$Estimate, 4)
  results_df$Lower_95_CI <- round(results_df$Lower_95_CI, 4)
  results_df$Upper_95_CI <- round(results_df$Upper_95_CI, 4)
  
  rownames(results_df) <- NULL
  
  print("--- Analysis Complete ---")
  return(results_df)
}

# #------------------------------------------------------------------------------#
# # EXECUTION EXAMPLE
# #------------------------------------------------------------------------------#
# data = data.sim(10) %>%
#   dplyr::rename(Time = Y, Event = D, Treatment = A) %>%
#   dplyr::select(-id)
# 
# summary(data)
# 
# df_results <- get_all_ATE_results(data, time.interest=10)
# 
# print("--------------------------------------------------")
# print(" FINAL ESTIMATED ATE RESULTS (95% CI) ")
# print("--------------------------------------------------")
# print(df_results)