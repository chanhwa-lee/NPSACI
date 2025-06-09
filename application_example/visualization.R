library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyr)
library(reshape2)
conflicts_prefer(dplyr::filter)

setwd("~/research/NPSACI/application_example")

###------------------------------------------------------###
###---------------- Estimation result  ------------------###
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
median_results_TypeB <- load_results("estimate/TypeB_nmax1000_r100_Scenario1/Rdata")
median_results_TPB   <- load_results("estimate/TPB_nmax1000_r100_Scenario1/Rdata")

# --- Risk over time plot ---
make_risk_plot = function(data, policy){
  
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
    labs(x = "Time (Days)", y = "Risk of Cholera", title = policy) +
    # scale_y_continuous(labels = function(y) y * 1000, limits = c(-0.0100,0.2185), breaks = 0.002*0:5) +
  
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

plot_risk_TypeB    <- make_risk_plot(median_results_TypeB, "TypeB")
plot_risk_TPB      <- make_risk_plot(median_results_TPB, "TPB")

plot_grid(plot_risk_TypeB, plot_risk_TPB, nrow = 1)


# --- Effects over theta plot ---
make_effect_plot = function(data, policy, times){
  
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

    # scale_y_continuous(labels = function(y) y*1000, limits = c(-0.0040,0.0015), breaks = c(-5:2)/1000) +

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
         title = policy) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.background = element_blank(),
          legend.key = element_blank())

}

# --- Make effect plots ---
plot_effect_TypeB    <- make_effect_plot(median_results_TypeB, "TypeB", c(180,360))
plot_effect_TPB      <- make_effect_plot(median_results_TPB, "TPB", c(180,360))

# --- Combined effect plots ---
plot_grid(plot_effect_TypeB + theme(legend.position = "none"),
          plot_effect_TPB,
          nrow = 2,
          rel_heights = c(8,9))
