###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)
library(reshape)
library(ggh4x)
library(scales)
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

data.plot_Ndist = data.frame()

for(ndist in 1:5){
  
  print(glue("Ndist: {ndist}"))
  
  # --- Load estimands and aggregate ---
  estimands <- load_results(glue("estimand/Ndist{ndist}")) %>% 
    dplyr::group_by(theta, tau, estimand) %>%
    dplyr::summarise(truth = mean(truth), .groups = "drop")
  
  # --- Load estimates and combine with estimands to compute metrics ---
  estimates <- load_results(glue("estimate/TypeB_m200_r100_Ndist{ndist}/Rdata"))
  
  result <- compute_metrics(estimates, estimands) %>% mutate(Ndist = ndist)
  
  data.plot_Ndist = rbind(data.plot_Ndist, result)
  
}


# Long form to generate plotting dataset
data.plot_Ndist = 
  reshape::melt(as.data.frame(data.plot_Ndist), 
                id.vars = c("estimand", "theta", "tau", "Ndist"), 
                var = "metric") %>%
  dplyr::filter(metric %in% c("Bias", "ESE", "COV", "UCOV"), 
                estimand %in% c("mu", "mu_1", "mu_0")) %>%
  dplyr::mutate(Ndist = factor(Ndist))

label_estimands = c("mu"   = expression(mu),
                    "mu_1" = expression(mu[1]),
                    "mu_0" = expression(mu[0]))

ggplot(data.plot_Ndist %>% 
         dplyr::filter(theta %in% c(0.3,0.4,0.5,0.6,0.7),
                       tau   %in% c(0.2,0.4)) %>%
         dplyr::mutate(theta = paste0( "alpha == ", theta),
                       tau   = glue("tau == {tau}")), 
       aes(x = estimand, y = value, group = Ndist, fill = Ndist))+ 
  
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
  
  scale_fill_discrete(name = "Distribution of N",
                      labels = c("1"  = TeX("N $\\in$\\{3,$\\cdots$,5\\}"),
                                 "2" = TeX("N $\\in$\\{5,$\\cdots$,20\\}"),
                                 "3" = TeX("N $\\in$\\{20,$\\cdots$,50\\}"),
                                 "4" = TeX("N $\\in$\\{50,$\\cdots$,100\\}"),
                                 "5" = TeX("N $\\sim$ Cholera Vaccine Data")
                                 )) +
  
  scale_x_discrete(labels = label_estimands) + 
  labs(x = "Estimand") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") + 
  ggh4x::facet_nested(metric ~ theta, scales = "free", labeller = label_parsed)

ggsave(filename = glue("FigS5.SimulationOver_Ndist.TypeB.pdf"), width = 9, height = 5)
