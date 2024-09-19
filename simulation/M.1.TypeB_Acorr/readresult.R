###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(reshape)

###------ Read estimand ------###
file.list <- list.files("estimand", pattern = "estimand.*rds", full.names = T)

estimands.list <- lapply(file.list, function(file) readRDS(file))

print(paste0(length(estimands.list), " estimand Rdata files were loaded"))

estimands = dplyr::bind_rows(estimands.list) %>% 
  dplyr::group_by(theta, tau, estimand) %>%
  dplyr::summarise(truth = mean(truth)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(theta = as.numeric(theta),
                tau = as.numeric(tau)) %>%
  dplyr::arrange(theta,
                 factor(estimand, levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")), 
                 tau)

print(estimands)

###------ Read estimator ------###

m = 200
r = 100

### Nonparametric estimation result ###

estimate.list.nonpara = list.files(paste0("estimate/m",m,"_r",r,"_nonpara/Rdata"), 
                           pattern = "estimate.*rds", full.names = T)

result.list.nonpara <- lapply(estimate.list.nonpara, function(file) readRDS(file)$result)

print(paste0(length(result.list.nonpara), " estimate Rdata files were loaded"))

## Aggregate estimation results
result.nonpara = bind_rows(result.list.nonpara, .id = "d")

## Merge with true estimands values
result.nonpara = result.nonpara %>%
  left_join(estimands, by = c("estimand", "theta", "tau")) %>%
  mutate(bias = est - truth,
         cov = I(abs(est - truth) < 1.96*se))

## Summarize performance metrics
result.table.nonpara = result.nonpara %>%
  dplyr::group_by(estimand, theta, tau) %>%
  dplyr::summarise(Truth = mean(truth, na.rm = T),
                   Bias = mean(bias, na.rm = T),
                   RMSE = sqrt(mean(bias^2, na.rm = T)),
                   ASE = mean(se, na.rm = T),
                   ESE = sd(est, na.rm = T),
                   Cov = round(mean(cov, na.rm = T),2)*100) %>%
  dplyr::filter(estimand != "te") %>%
  dplyr::arrange(theta,
                 factor(estimand, levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")), 
                 tau)

result_combined = 
inner_join(result.table.nonpara %>% filter(tau == 0.3, theta %in% c(0.3,0.5,0.7)),
           result.table.nonpara %>% filter(tau == 0.5, theta %in% c(0.3,0.5,0.7)),
           by = c("estimand", "theta"),
           suffix = c(".1", ".2"))

# Identify continuous columns (numeric or integer)
continuous_cols <- sapply(result_combined, is.numeric)

# Mutiplying by 100
result_combined_multiplied = result_combined
multiplying_cols <- c("Truth.1", "Bias.1", "RMSE.1", "ASE.1", "ESE.1", "Truth.2", "Bias.2", "RMSE.2", "ASE.2", "ESE.2")
result_combined_multiplied[multiplying_cols] <- round(result_combined_multiplied[multiplying_cols]*100, 1)

result_combined_multiplied %>% print(n=100)

write.csv(result_combined_multiplied,
          paste0("Table3.simulation_TypeB_Acorr_m",m,"_r",r,"_taus_multiplied.csv"))


# Round up to 3 decimals for continuous columns
result_combined[continuous_cols] <- round(result_combined[continuous_cols], 3)

result_combined %>% print(n=100)

write.csv(result_combined, 
          paste0("Table3.simulation_TypeB_Acorr_m",m,"_r",r,"_taus.csv"))
