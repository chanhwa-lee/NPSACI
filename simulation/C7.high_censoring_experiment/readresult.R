###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(reshape)
library(ggplot2)
library(conflicted)
library(glue)
library(cowplot)
library(latex2exp)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::mutate)

# --- Function to load and process result files ---
load_results <- function(path) {
  
  file_list <- list.files(path, pattern = ".*rds", full.names = TRUE)
  results <- lapply(file_list, function(f) readRDS(f))
  
  message(length(file_list), " Rdata files were loaded from ", path)
  if(length(file_list) < 10) return()
  
  results <- bind_rows(results, .id = "d")
  
  results <- results %>%
    mutate(estimand = factor(estimand, 
                             levels = c("mu", "mu_1", "mu_0", "de", "oe", "se_1", "se_0"),
                             ordered = TRUE),
           theta = round(theta, 4),
           tau = round(tau, 4))
  
  return(results)
}

# --- Function to compute performance metric ---
compute_metrics <- function(estimates, estimands){
  
  estimates = estimates %>%
    filter(abs(est) <= 1) %>%  ## Filter out extreme estimates
    left_join(estimands, by = c("estimand", "theta", "tau")) %>%
    mutate(bias = est - truth,
           cov  = I(PCL <= truth & truth <= PCU),
           ucov = I(UCL <= truth & truth <= UCU))
  
  result_UCov = estimates %>%
    dplyr::group_by(estimand, tau, d) %>%
    dplyr::summarize(ucov = mean(ucov, na.rm = TRUE) >= 1, .groups = "drop_last") %>%
    dplyr::summarize(UCOV = round(mean(ucov, na.rm = T)*1000)/10, .groups = "drop")
  
  result = estimates %>%
    dplyr::group_by(estimand, theta, tau) %>%
    dplyr::summarise(Truth = mean(truth, na.rm = T),
                     Bias = mean(bias, na.rm = T),
                     AbsBias = mean(abs(bias), na.rm = T),
                     RelBias = abs(mean(bias, na.rm = T) / Truth),
                     RMSE = sqrt(mean(bias^2, na.rm = T)),
                     ASE = mean(se, na.rm = T),
                     ESE = sd(est, na.rm = T),
                     COV = round(mean(cov, na.rm = T)*1000)/10,
                     .groups = "drop") %>%
    left_join(result_UCov, by = c("estimand", "tau"))
  
  return(result)
  
}

# --- Load estimands and estimates ---
estimands <- load_results("estimand/TypeB") %>% 
  dplyr::group_by(theta, tau, estimand) %>%
  dplyr::summarise(truth = mean(truth), .groups = "drop")

df_result = data.frame()

for(censoring_rate in c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)){
  
  print(glue("Censoring rate: {censoring_rate}"))
  
  estimates <- load_results(glue("estimate/TypeB_m200_r100_CensoringRate{censoring_rate}/Rdata"))
  
  if(is.null(estimates)) {
    print("Estimates not loaded")
    next
  }
  
  result <- compute_metrics(estimates, estimands) %>%
    filter(tau %in% c(0.2,0.4),
           estimand %in% c("mu", "mu_1", "mu_0"),
           theta %in% c(0.3, 0.45, 0.6) ) %>%
    select(estimand, theta, tau, Truth, Bias, RMSE, ASE, ESE, COV, UCOV) %>%
    arrange(tau, theta) %>%
    mutate(censoring_rate = censoring_rate) %>% 
    print(n=100)
  
  df_result = rbind(df_result, result)
  
  write.csv(table, glue("Table.Simulation_TypeB_CensoringRate{censoring_rate}.csv"))
}

# Long form to generate plotting dataset
data.plot_censoring = 
  reshape::melt(as.data.frame(df_result), 
                id.vars = c("estimand", "theta", "tau", "censoring_rate"), 
                var = "metric") %>%
  dplyr::filter(metric %in% c("Bias", "ESE", "COV", "UCOV"), 
                estimand %in% c("mu", "mu_1", "mu_0")) %>%
  dplyr::mutate(censoring_rate = factor(censoring_rate))

label_estimands = c("mu"   = expression(mu),
                    "mu_1" = expression(mu[1]),
                    "mu_0" = expression(mu[0]))

ggplot(data.plot_censoring %>% 
         dplyr::filter(theta %in% c(0.3,0.45,0.6),
                       tau   %in% c(0.2,0.4)) %>%
         dplyr::mutate(theta = paste0( "alpha == ", theta),
                       tau   = glue("tau == {tau}")), 
       aes(x = estimand, y = value, group = censoring_rate, fill = censoring_rate))+ 
  
  geom_bar(stat = "identity", position = position_dodge()) +
  
  ### 95% CI benchmark line ###
  geom_hline(data = expand.grid(metric = "COV"), 
             aes(yintercept = 95),
             color = "black",
             linetype = "dashed") +
  
  ### 95% UCB benchmark line ###
  geom_hline(data = expand.grid(metric = "UCOV"), 
             aes(yintercept = 95),
             color = "black",
             linetype = "dashed") +
  
  scale_fill_discrete(name = "Censoring Rate", guide = guide_legend(nrow = 1)) +
  
  scale_x_discrete(labels = label_estimands) + 
  labs(x = "Estimand") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") + 
  ggh4x::facet_nested(metric ~ theta + tau, scales = "free", labeller = label_parsed)

ggsave(filename = glue("FigS.SimulationOver_Censoring_Rate.TypeB.pdf"), width = 9, height = 5)
