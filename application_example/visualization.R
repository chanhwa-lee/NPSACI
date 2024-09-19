library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(reshape2)

setwd("~/research/NPSACI/application_example")

###------------------------------------------------------###
###---------------- Estimation result  ------------------###
###------------------------------------------------------###

###--- TypeB results ---###

## Read estimator

result.list.TypeB <- list.files(paste0("TypeB/estimate/S5_r100_nmax10_nonpara/Rdata"), 
                          pattern = "estimate.*rds", full.names = T)

M <- length(result.list.TypeB)

result.TypeB <- list()

for(i in 1:M){
  result.TypeB[[i]] = readRDS(result.list.TypeB[i])$result
}

print(paste0(M, " estimate Rdata files were loaded"))

## Aggregate estimation results

result.TypeB = bind_rows(result.TypeB, .id = "s")

result.median.TypeB = result.TypeB %>%
  group_by(estimand, theta, tau) %>%
  summarise(est = median(est, na.rm = T),
            se  = median(se, na.rm = T)) %>%
  ungroup() %>%
  mutate(estimand = factor(estimand, 
                           levels = c("mu", "mu_1", "mu_0", "de", "oe", "se_1", "se_0", "te"), 
                           ordered = TRUE)) %>%
  arrange(theta, estimand)

###--- TPB results ---###

## Read estimator

result.list.TPB <- list.files(paste0("TPB/estimate/S5_r100_nmax10_nonpara/Rdata"), 
                              pattern = "estimate.*rds", full.names = T)

M <- length(result.list.TPB)

result.TPB <- list()

for(i in 1:M){
  result.TPB[[i]] = readRDS(result.list.TPB[i])$result
}

print(paste0(M, " estimate Rdata files were loaded"))


## Aggregate estimation results

result.TPB = bind_rows(result.TPB, .id = "s")

result.median.TPB = result.TPB %>%
  group_by(estimand, theta, tau) %>%
  summarise(est = median(est, na.rm = T),
            se  = median(se, na.rm = T)) %>%
  ungroup() %>%
  mutate(estimand = factor(estimand, 
                           levels = c("mu", "mu_1", "mu_0", "de", "oe", "se_1", "se_0", "te"), 
                           ordered = TRUE)) %>%
  arrange(theta, estimand)



###------ Risk over time plot ------###

## Type B ##

p.risk.TypeB = 
  ggplot(data = result.median.TypeB %>% 
           filter(estimand %in% c("mu_1", "mu_0")) %>%
           filter(tau <= 500) %>%
           filter(theta %in% c(0.3, 0.45, 0.6)) %>%
           mutate(theta = paste0("alpha == ", theta)),
         aes(x = tau, y = est, group = estimand, color = estimand)) +
  geom_line(aes(linetype = estimand)) + 
  facet_grid(. ~ theta, scales = "free_y", labeller = label_parsed) +
  
  labs(title = "Type B") + 
  
  # Change x-axis and y-axis labels
  labs(x = "Time (Days)", y = "Risk of Cholera") +
  
  scale_y_continuous(limits = c(0,0.26)) +
  
  # Change line type for 'mu_1' and 'mu_0' and customize legend labels
  scale_linetype_manual(values = c("mu_0" = "solid", "mu_1" = "dashed"), 
                        labels = c("mu_0" = "Unvaccinated", "mu_1" = "Vaccinated"),
                        guide = guide_legend(title = NULL)) +
  
  # Change color for 'mu_1' and 'mu_0' and customize legend labels
  scale_color_manual(values = c("mu_0" = "black", "mu_1" = "blue"), 
                     labels = c("mu_0" = "Unvaccinated", "mu_1" = "Vaccinated"),
                     guide = guide_legend(title = NULL)) +
  
  theme_bw()+ 
  
  theme(legend.position="bottom")

p.risk.TypeB

## TPB ##

p.risk.TPB = 
  ggplot(data = result.median.TPB %>% 
           filter(estimand %in% c("mu_1", "mu_0")) %>%
           filter(tau <= 500) %>%
           filter(theta %in% c(0, 0.25, 0.5)) %>%
           mutate(theta = paste0("rho == ", theta)),
         aes(x = tau, y = est, group = estimand, color = estimand)) +
  geom_line(aes(linetype = estimand)) + 
  
  facet_grid(. ~ theta, labeller = label_parsed) +
  
  labs(title = "TPB") + 
  
  # Change x-axis and y-axis labels
  labs(x = "Time (Days)", y = "Risk of Cholera") +
  
  scale_y_continuous(limits = c(0,0.26)) +
  
  # Change line type for 'mu_1' and 'mu_0' and customize legend labels
  scale_linetype_manual(values = c("mu_0" = "solid", "mu_1" = "dashed"), 
                        labels = c("mu_0" = "Unvaccinated", "mu_1" = "Vaccinated"),
                        guide = guide_legend(title = NULL)) +
  
  # Change color for 'mu_1' and 'mu_0' and customize legend labels
  scale_color_manual(values = c("mu_0" = "black", "mu_1" = "blue"), 
                     labels = c("mu_0" = "Unvaccinated", "mu_1" = "Vaccinated"),
                     guide = guide_legend(title = NULL)) +
  
  theme_bw() + 
  theme(legend.position="bottom") 


p.risk.TPB

## Risk plot combined (TypeB, TPB)
plot_grid(
  plot_grid(p.risk.TypeB + theme(legend.position="none"), 
            p.risk.TPB + theme(legend.position="none"), 
            nrow = 1,
            align = "v"),
  get_legend(
    p.risk.TypeB
  ),
  nrow = 2,
  rel_heights = c(2,0.1)
)



###------ Effects over thetas ------###

## Type B ##

p.effects.TypeB = 
  ggplot(data = result.median.TypeB %>% 
           dplyr::filter(estimand %in% c("de", "se_1", "se_0", "oe")) %>% 
           mutate(estimand = case_when(
             estimand == "de"   ~ 'DE["B"](tau:alpha)',
             estimand == "se_1" ~ 'SE["B,1"](tau:alpha,0.45)',
             estimand == "se_0" ~ 'SE["B,0"](tau:alpha,0.45)',
             estimand == "oe"   ~ 'OE["B"](tau:alpha,0.45)'
           )) %>% 
           filter(tau %in% c(250, 500)) %>%
           filter(0.3 <= theta, theta <= 0.6) %>%
           # filter(round(theta, digits = 3) %in% round(seq(0.3, 0.6, length.out = 61), digits = 3)) %>%
           mutate(tau = paste0("tau == ", tau)), 
         aes(x = theta, y = est, group = tau)) +
  geom_errorbar(aes(ymin = est - 1.96*se, ymax = est + 1.96*se), 
                size = 0.4, width = 0.00, col="#997570") +
  geom_line(col = "red") + 
  geom_point(size = 0.1, col = "red") +
  geom_hline(yintercept = 0, color = "blue", size = 0.5, linetype = "dashed") + 
  scale_y_continuous(labels = function(y) y*1000) + 
  facet_grid(tau ~ estimand, labeller = label_parsed) + 
  # Change x-axis and y-axis labels
  # labs(title = "Estimated Causal Effects under Type B policy")+
  labs(title = "Type B") + 
  labs(x = expression(alpha), y = "Effects") +
  theme_bw()

p.effects.TypeB


## TPB ##

p.effects.TPB = 
  ggplot(data = result.median.TPB %>% 
           dplyr::filter(estimand %in% c("de", "se_1", "se_0", "oe")) %>% 
           mutate(estimand = case_when(
             estimand == "de"   ~ 'DE["TPB"](tau:rho)',
             estimand == "se_1" ~ 'SE["TPB,1"](tau:rho,0)',
             estimand == "se_0" ~ 'SE["TPB,0"](tau:rho,0)',
             estimand == "oe"   ~ 'OE["TPB"](tau:rho,0)'
           )) %>% 
           filter(tau %in% c(250, 500)) %>%
           filter(theta <= 0.5) %>%
           mutate(tau = paste0("tau == ", tau)), 
         aes(x = theta, y = est, group = tau)) +
  geom_errorbar(aes(ymin = est - 1.96*se, ymax = est + 1.96*se), 
                size = 0.4, width = 0.00, col="#997575") +
  geom_line(col = "red") + 
  geom_point(size = 0.1, col = "red") +
  geom_hline(yintercept = 0, color = "blue", size = 0.5, linetype = "dashed") + 
  scale_y_continuous(labels = function(y) y*1000) + 
  facet_grid(tau ~ estimand, labeller = label_parsed) + 
  labs(title = "TPB") + 
  labs(x = expression(rho), y = "Effects") +
  theme_bw()

p.effects.TPB

## Effects plot combined (TypeB, TPB)
plot_grid(p.effects.TypeB, p.effects.TPB, ncol = 1)
