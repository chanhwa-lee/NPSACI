library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyr)
library(reshape2)
conflicts_prefer(dplyr::filter)

setwd("~/research/NPSACI/application")

###------------------------------------------------------###
###---------------- Estimation result  ------------------###
###-------------- (Main text Section 6) -----------------###
###------------------------------------------------------###


# --- Function to load and process result files ---
load_results <- function(path) {
  file_list <- list.files(path, pattern = "estimate.*rds", full.names = TRUE)
  
  results <- lapply(file_list, function(f) readRDS(f)$result)
  
  message(length(file_list), " estimate Rdata files were loaded from ", path)
  
  results <- bind_rows(results, .id = "s")
  
  median_results <- results %>%
    group_by(estimand, theta, tau) %>%
    summarise(across(c(est, se, PCL, PCU, UCL, UCU), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(
      estimand = factor(estimand, 
                        levels = c("mu", "mu_1", "mu_0", "de", "oe", "se_1", "se_0", "te"),
                        ordered = TRUE),
      theta = round(theta, 4)
    ) %>%
    arrange(theta, estimand)
  
  return(median_results)
}

# --- Load median results for each combination ---
median_results_TypeB_nonpara <- load_results("estimate/TypeB_nmax1000_r100_Scenario1/Rdata")
median_results_TypeB_para    <- load_results("estimate/TypeB_nmax1000_r100_Scenario43/Rdata")
median_results_TPB_nonpara   <- load_results("estimate/TPB_nmax1000_r100_Scenario1/Rdata")
median_results_TPB_para      <- load_results("estimate/TPB_nmax1000_r100_Scenario43/Rdata")

# --- Risk over time plot ---
make_risk_plot = function(data, policy, method){
  
  # Set theta values depending on policy
  if (policy == "TypeB") {
    thetas <- c(0.3, 0.45, 0.6)
  } else if (policy == "TPB") {
    thetas <- c(0, 0.25, 0.5)
  } else {
    stop("Unknown policy type")
  }
  
  p = ggplot(data = data %>%
          filter(estimand %in% c("mu_1", "mu_0")) %>%
          filter(tau <= 460) %>%
          filter(theta %in% thetas) %>%
          mutate(theta = paste0( ifelse(policy == "TypeB", "alpha == ", "rho == "), theta)),
         aes(x = tau, y = est, group = estimand, color = estimand)) +
    geom_line(aes(linetype = estimand)) +
    facet_grid(. ~ theta, labeller = label_parsed) +
    labs(x = "Time (Days)", y = "Risk of Cholera", title = paste(policy, "-", method)) +
    scale_y_continuous(labels = function(y) y * 1000, limits = c(0,0.0085), breaks = 0.002*0:5) +
  
    # Change line type for 'mu_1' and 'mu_0' and customize legend labels
    scale_linetype_manual(values = c("mu_0" = "solid", "mu_1" = "dashed"),
                          labels = c("mu_0" = "Unvaccinated", "mu_1" = "Vaccinated"),
                          guide = guide_legend(title = NULL)) +

    # Change color for 'mu_1' and 'mu_0' and customize legend labels
    scale_color_manual(values = c("mu_0" = "black", "mu_1" = "blue"),
                       labels = c("mu_0" = "Unvaccinated", "mu_1" = "Vaccinated"),
                       guide = guide_legend(title = NULL)) +
    theme_bw() +
    theme(legend.position = "bottom")

  return(p)
}


# --- Make risk plots ---
plot_risk_TypeB_nonpara    <- make_risk_plot(median_results_TypeB_nonpara, "TypeB", "Nonparametric")
plot_risk_TypeB_para       <- make_risk_plot(median_results_TypeB_para, "TypeB", "Parametric")
plot_risk_TPB_nonpara      <- make_risk_plot(median_results_TPB_nonpara, "TPB", "Nonparametric")
plot_risk_TPB_para         <- make_risk_plot(median_results_TPB_para, "TPB", "Parametric")


# --- Combined risk plots ---
plot_grid(plot_risk_TypeB_nonpara,
          plot_risk_TPB_nonpara  ,
          plot_risk_TypeB_para   ,
          plot_risk_TPB_para     ,
          nrow = 2)

plot_grid(plot_risk_TypeB_nonpara + labs(title = "TypeB"),
          plot_risk_TPB_nonpara + labs(title = "TPB"),
          nrow = 1)

ggsave("Fig2.Risk_over_time_combined.pdf", width = 9, height = 3.5)

# --- Effects over theta plot ---
make_effect_plot = function(data, policy, method, times){
  
  # Set theta values depending on policy
  if (policy == "TypeB") {
    thetas <- c(0.3, 0.45, 0.6)
    ## Define label expressions using a named list
    label_greek <- as_labeller(c(
      "mu"   = "mu['B'](alpha)", 
      "mu_1" = "mu['B,1'](alpha)", 
      "mu_0" = "mu['B,0'](alpha)", 
      "de"   = 'DE["B"](tau:alpha)',
      "se_1" = 'SE["B,1"](tau:alpha,0.45)',
      "se_0" = 'SE["B,0"](tau:alpha,0.45)',
      "oe"   = 'OE["B"](tau:alpha,0.45)'
    ), label_parsed)
    
  } else if (policy == "TPB") {
    thetas <- c(0, 0.25, 0.5)
    ## Define label expressions using a named list
    label_greek <- as_labeller(c(
      "mu"   = "mu['TPB'](rho)", 
      "mu_1" = "mu['TPB,1'](rho)", 
      "mu_0" = "mu['TPB,0'](rho)", 
      "de"   = 'DE["TPB"](tau:rho)',
      "se_1" = 'SE["TPB,1"](tau:rho,0)',
      "se_0" = 'SE["TPB,0"](tau:rho,0)',
      "oe"   = 'OE["TPB"](tau:rho,0)'
    ), label_parsed)
    
  } else {
    stop("Unknown policy type")
  }
  
  ggplot(data %>%
           dplyr::filter(estimand %in% c("de", "se_1", "se_0", "oe")) %>%
           dplyr::mutate(estimand = factor(estimand, levels = c("de", "oe", "se_0", "se_1"))) %>%
           dplyr::filter(tau %in% times) %>%
           mutate(tau = factor(paste0("tau == ", tau), levels = paste0("tau == ", times))),
         aes(x = theta, y = est)) +
    ## Uniform confidence band
    geom_ribbon(aes(ymin = UCL, ymax = UCU, fill = "95% Uniform Confidence Band"), alpha = 1) +

    ## Pointwise confidence interval
    geom_ribbon(aes(ymin = PCL, ymax = PCU, fill = "95% Pointwise Confidence Interval"), alpha = 1) +

    ## Estimate line
    geom_line(aes(color = "Estimate"), linewidth = 0.5) +

    ### Effect = 0 benchmark line ###
    geom_hline(data = expand.grid(estimand = c("de", "se_1", "se_0", "oe")),
               aes(yintercept = 0),
               color = "red",
               linetype = "dashed") +

    ## Facet
    facet_grid(tau ~ estimand, labeller = labeller(
      estimand = label_greek,
      tau = label_parsed
    )) +

    scale_y_continuous(labels = function(y) y*1000, limits = c(-0.0040,0.0015), breaks = c(-5:2)/1000) +

    ## Manual legend control
    scale_fill_manual(name = NULL,
                      values = c("95% Pointwise Confidence Interval" = "lightblue3",
                                 "95% Uniform Confidence Band" = "lightblue1"),
                      guide = guide_legend(title = NULL)) +
    scale_color_manual(name = NULL,
                       values = c("Estimate" = "black"),
                       guide = guide_legend(title = NULL)) +

    ## Labels and theme
    labs(x = ifelse(policy == "TypeB", expression(alpha), expression(rho)),
         y = "Effects",
         title = paste(policy, "-", method)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.background = element_blank(),
          legend.key = element_blank())

}

# --- Make effect plots ---
plot_effect_TypeB_nonpara    <- make_effect_plot(median_results_TypeB_nonpara, "TypeB", "Nonparametric", c(180,360))
plot_effect_TypeB_para       <- make_effect_plot(median_results_TypeB_para, "TypeB", "Parametric", c(180,360))
plot_effect_TPB_nonpara      <- make_effect_plot(median_results_TPB_nonpara, "TPB", "Nonparametric", c(180,360))
plot_effect_TPB_para         <- make_effect_plot(median_results_TPB_para, "TPB", "Parametric", c(180,360))

# --- Combined effect plots ---
plot_grid(plot_effect_TypeB_nonpara,
          plot_effect_TPB_nonpara  ,
          plot_effect_TypeB_para   ,
          plot_effect_TPB_para     ,
          nrow = 2)

plot_grid(plot_effect_TypeB_nonpara + labs(title = "TypeB") + theme(legend.position = "none"),
          plot_effect_TPB_nonpara   + labs(title = "TPB"),
          nrow = 2,
          rel_heights = c(8,9))

ggsave("Fig3.Effects_over_theta_combined.pdf", width = 9, height = 8)




###--- Finder function ---###

query = function(data, estimand_, theta_, tau_){
  
  data %>% 
    filter(estimand == estimand_, theta == theta_, tau == tau_) %>% 
    mutate(est = est*1000, se = se*1000) %>%
    mutate(lo = round(est - 1.96*se,4),
           up = round(est + 1.96*se,4)) %>%
    as.data.frame()
  
}

query(data = median_results_TypeB_nonpara, estimand_ = "de", theta_ = 0.3, tau_ = 360)

query(data = median_results_TPB_nonpara, estimand_ = "se_0", theta_ = 0.5, tau_ = 360)

query(data = median_results_TypeB_nonpara, estimand_ = "oe", theta_ = 0.6, tau_ = 360)




###------------------------------------------------------###
###---- Effects over thetas - multiple time points ------###
###----------- (Supplementary Section D.2) --------------###
###------------------------------------------------------###

## Type B ##
make_effect_plot(median_results_TypeB_nonpara, "TypeB", "Nonparametric", c(90,180,270,360,450)) +
  scale_y_continuous(labels = function(y) y*1000, limits = c(-0.0065,0.0030), breaks = c(-6,-4,-2,0,2)/1000) + 
  labs(title = "TypeB")

ggsave("FigS9.Effects_over_theta_TypeB_multiple_timepoints.pdf", width = 8, height = 7)

## TPB ##
make_effect_plot(median_results_TPB_nonpara, "TPB", "Nonparametric", c(90,180,270,360,450)) +
  scale_y_continuous(labels = function(y) y*1000, limits = c(-0.0060,0.0015), breaks = c(-6,-4,-2,0,2)/1000) + 
  labs(title = "TPB")

ggsave("FigS10.Effects_over_theta_TPB_multiple_timepoints.pdf", width = 8, height = 7)





###------------------------------------------------------###
###------------------- Choice of S ----------------------###
###----------- (Supplementary Section D.3) --------------###
###------------------------------------------------------###

S.values = c(1,3,5,10,15)

file_list <- list.files("./estimate/TypeB_nmax1000_r100_Scenario1/Rdata", pattern = "estimate.*rds", full.names = TRUE)

result.TypeB <- lapply(file_list, function(f) readRDS(f)$result)
result.TypeB <- bind_rows(result.TypeB, .id = "s")

result.median.TypeB.S = data.frame()

for(S in S.values){
  result.median.TypeB.temp = result.TypeB %>%
    mutate(s = as.numeric(s)) %>%
    filter(s <= S, estimand == "mu") %>%
    group_by(theta, tau) %>%
    summarise(est = median(est, na.rm = T), .groups = "drop") %>%
    arrange(theta, tau) %>%
    mutate("S" = S)
  
  result.median.TypeB.S = rbind(result.median.TypeB.S, result.median.TypeB.temp)
  
}

result.median.TypeB.S$S = as.factor(result.median.TypeB.S$S)

ggplot(data = result.median.TypeB.S %>% 
           filter(tau <= 470) %>%
           filter(round(theta,4) %in% c(0.3, 0.45, 0.6)) %>%
           mutate(theta = paste0("alpha == ", theta)),
         aes(x = tau, y = est, group = S, color = S)) +
  geom_line(aes(linetype = S)) + 
  facet_grid(. ~ theta, scales = "free_y", labeller = label_parsed) +
  labs(x = "Time (Days)", y = "Risk of Cholera") +
  scale_y_continuous(limits = c(0,0.008), labels = function(y) y*1000) +
  theme_bw()

ggsave("FigS11.ScompTypeB.pdf", width = 8, height = 3)
