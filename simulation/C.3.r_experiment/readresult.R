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
  
  ## Summarize performance metrics
  result_table = estimates %>%
    # dplyr::filter(estimand %in% c("mu_1", "mu_0")) %>%
    dplyr::group_by(estimand, theta, tau) %>%
    dplyr::summarise(Truth = mean(truth, na.rm = T),
                     Bias = mean(bias, na.rm = T),
                     RMSE = sqrt(mean(bias^2, na.rm = T)),
                     ASE = mean(se, na.rm = T),
                     ESE = sd(est, na.rm = T),
                     Cov = mean(cov, na.rm = T)) %>%
    dplyr::arrange(theta,
                   factor(estimand, levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")), 
                   tau)
  
  return(result_table)
}

## Read estimates over r values ##

data.plot = data.frame()

for(r in c(5,10,20,50,100,200,500)){
  
  print(paste0("r = ",r))
  
  estimate.list = list.files(paste0("estimate/m100_r",r,"_nonpara/Rdata"), 
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
    mutate(r = r)
  
  data.plot = rbind(data.plot, result_table)
  
}

## Visualization ##

library(tidyverse)

data.plot.2 <- data.plot %>%
  mutate(SER = ASE / ESE) %>%
  mutate(Cov = round(Cov, 2)*100) %>%
  pivot_longer(cols = c(Bias, RMSE, ASE, ESE, Cov, SER),
               names_to = "type",
               values_to = "value")

### Bar chart ###

ggplot(data = data.plot.2 %>% 
         filter(theta >= 0.3, theta <= 0.7) %>%
         filter(type %in% c("Bias", "ESE", "Cov")) %>%
         mutate(type = factor(type, levels = c("Bias", "ESE", "Cov"), ordered = T)) %>%
         mutate(delta = paste0("alpha == ", theta)),
       aes(x = estimand, y = value, group = r, fill = r)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ### 95% CI benchmark line ###
  geom_hline(data = expand.grid(r = unique(data.plot.2$r),
                                type = "Cov",
                                theta = paste0("alpha == ", 0.1*(3:7))), 
             aes(yintercept = 95),
             color = "blue",
             linetype = "dashed") +
  facet_grid(type ~ delta, scales = "free", labeller = label_parsed) +
  scale_x_discrete(labels = c("mu" = expression(mu(alpha)),
                              "mu_1" = expression(mu[1](alpha)),
                              "mu_0" = expression(mu[0](alpha)))) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

ggsave(filename = "FigS2.r.comparison.m100.D1000.pdf", width = 8, height = 4)
