###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(latex2exp)
library(tidyverse)

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
  
  # # Identify continuous columns (numeric or integer)
  # continuous_cols <- sapply(result_table, is.numeric)
  # 
  # # Round up to 3 decimals for continuous columns
  # result_table[continuous_cols] <- round(result_table[continuous_cols], 3)
  
  return(result_table)
}



data.plot = data.frame()

cases = c("N3",
          "N3_5",
          "N5",
          "N5_20",
          "N20",
          "N20_50",
          "N50_100")

for(case in cases){
  
  ###------ Read estimand ------###
  file.list <- list.files(paste0(case,"/estimand"), pattern = "estimand.*rds", full.names = T)
  
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
  
  estimate.list = list.files(paste0(case,"/estimate/m100_r100_nonpara/Rdata"), 
                             pattern = "estimate.*rds", full.names = T)
  
  estimates <- lapply(estimate.list, function(file) readRDS(file)$result)
  
  print(paste0(length(estimates), " estimate Rdata files were loaded"))
  
  ## Aggregate estimation results
  estimates = bind_rows(estimates, .id = "d")
  
  result_table = 
    result.summarise(estimates, estimands) %>% 
    filter(estimand %in% c("mu", "mu_1", "mu_0")) %>%
    mutate(case = case)
  
  data.plot = rbind(data.plot, result_table)
  
}




## Visualization ##
data.plot.2 <- data.plot %>%
  mutate(SER = ASE / ESE) %>%
  mutate(Cov = round(Cov, 2)*100) %>%
  pivot_longer(cols = c(Bias, RMSE, ASE, ESE, Cov, SER),
               names_to = "type",
               values_to = "value")

### Bar chart ###
ggplot(data = data.plot.2 %>% 
         filter(tau == 0.3) %>%
         filter(theta >= 0.3, theta <= 0.7) %>% 
         filter(type %in% c("Bias", "ESE", "Cov")) %>%
         mutate(case = factor(case, levels = cases, ordered = T)) %>%
         mutate(type = factor(type, levels = c("Bias", "ESE", "Cov"), ordered = T)) %>%
         mutate(theta = paste0("alpha == ", theta)),
       aes(x = estimand, y = value, group = case, fill = case)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ### 95% CI benchmark line ###
  geom_hline(data = expand.grid(method = cases,
                                type = "Cov",
                                theta = paste0("alpha == ", 0.1*(3:7))), 
             aes(yintercept = 95),
             color = "blue",
             linetype = "dashed") +
  facet_grid(type ~ theta, scales = "free", labeller = label_parsed) +
  scale_x_discrete(labels = c("mu" = expression(mu(alpha)),
                              "mu_1" = expression(mu[1](alpha)),
                              "mu_0" = expression(mu[0](alpha)))) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  scale_fill_discrete(name = "Distribution of N",
                      labels = c("N3"    = TeX("N $\\equiv$ 3"),
                                 "N3_5"  = TeX("N $\\in$\\{3,$\\cdots$,5\\}"),
                                 "N5"    = TeX("N $\\equiv$ 5"),
                                 "N5_20" = TeX("N $\\in$\\{5,$\\cdots$,20\\}"),
                                 "N20"    = TeX("N $\\equiv$ 20"),
                                 "N20_50" = TeX("N $\\in$\\{20,$\\cdots$,50\\}"),
                                 "N50_100" = TeX("N $\\in$\\{50,$\\cdots$,100\\}")))

ggsave(filename = "FigS4.Ndist.comparison.m100.D1000.pdf", width = 8, height = 5)
