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

## Help function ##
result.summarise = function(estimates, estimands){
  
  ## Merge with true estimands values
  estimates = estimates %>%
    left_join(estimands, by = c("estimand", "theta", "tau")) %>%
    mutate(bias = est - truth,
           cov = I(abs(est - truth) < 1.96*se))
  
  estimates = estimates %>% filter(abs(est) < 1)
  
  ## Summarize performance metrics
  result_table = estimates %>%
    # dplyr::filter(estimand %in% c("mu_1", "mu_0")) %>%
    dplyr::group_by(estimand, theta, tau) %>%
    dplyr::summarise(Truth = mean(truth, na.rm = T),
                     Bias = mean(bias, na.rm = T),
                     RMSE = sqrt(mean(bias^2, na.rm = T)),
                     ASE = mean(se, na.rm = T),
                     ESE = sd(est, na.rm = T),
                     Cov = round(mean(cov, na.rm = T),2)*100) %>%
    dplyr::arrange(theta,
                   factor(estimand, levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")), 
                   tau)
  
  # # Identify continuous columns (numeric or integer)
  # continuous_cols <- sapply(result_table, is.numeric)
  # 
  # # Round up to 3 decimals for continuous columns
  # result_table[continuous_cols] <- round(result_table[continuous_cols], 3)
  
  return(result_table)
}

## Read estimates over r values ##

data.plot = data.frame()

for(method in c("BDD","noBDD")){
  
  print(paste0("method = ", method))
  
  estimate.list = list.files(paste0("estimate/m200_r100_nonpara_",method,"/Rdata"), 
                             pattern = "estimate.*rds", full.names = T)
  
  estimates <- lapply(estimate.list, function(file) readRDS(file)$result)
  
  print(paste0(length(estimates), " estimate Rdata files were loaded"))
  
  ## Aggregate estimation results
  estimates = bind_rows(estimates, .id = "d")
  
  result_table = 
    result.summarise(estimates, estimands) %>% 
    filter(estimand %in% c("mu", "mu_1", "mu_0")) %>%
    # filter(theta %in% c(0.4,0.5,0.6)) %>%
    filter(tau == 0.3) %>%
    mutate(method = method)
  
  data.plot = rbind(data.plot, result_table)
  
}

## Visualization ##

library(tidyverse)

data.plot.2 <- data.plot %>%
  mutate(SER = ASE / ESE) %>%
  pivot_longer(cols = c(Bias, RMSE, ASE, ESE, Cov, SER),
               names_to = "type",
               values_to = "value")

### Bar chart ###
ggplot(data = data.plot.2 %>% 
         filter(theta >= 0.3, theta <= 0.7) %>%
         filter(type %in% c("Bias", "ESE", "Cov")) %>%
         mutate(type = factor(type, levels = c("Bias", "ESE", "Cov"), ordered = T)) %>%
         mutate(theta = paste0("alpha == ", theta)),
       aes(x = estimand, y = value, group = method, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(data = expand.grid(method = c("BDD", "noBDD"),
                                type = "Cov",
                                theta = paste0("alpha == ", 0.1*(3:7))), 
             aes(yintercept = 95),
             color = "blue",
             linetype = "dashed") +
  facet_grid(type ~ theta, scales = "free", labeller = label_parsed) +
  scale_x_discrete(labels = c("mu" = expression(mu(alpha)),
                              "mu_1" = expression(mu[1](alpha)),
                              "mu_0" = expression(mu[0](alpha)))) +
  scale_fill_discrete("Method", label = c("BDD" = "Bounded", "noBDD" = "Unbounded")) + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom") 

ggsave(filename = "FigS1.BDD.comparison.m200.D1000.pdf", width = 8, height = 4)
