###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.5")
library(dplyr)
library(reshape2)
library(ggh4x)
library(scales)
library(ggplot2)
library(conflicted)
library(glue)
library(cowplot)
library(latex2exp)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(glue)


conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::mutate)

m=10
N_values <- c(100, 300, 500, 700, 900) # Add your desired Ns here

# m=200
# N_values <- c(20, 50, 100, 200) # Add your desired Ns here

# --- Function to load and process result files ---
load_results <- function(path) {
  file_list <- list.files(path, pattern = ".*rds", full.names = TRUE)
  if(length(file_list) == 0) return(NULL)
  
  results <- lapply(file_list, function(f) readRDS(f))
  message(length(file_list), " files loaded from ", path)
  
  if(length(results) < 1) return(NULL)
  
  results <- bind_rows(results, .id = "d")
  results <- results %>%
    dplyr::filter(estimand %in% c("mu", "mu_1", "mu_0")) %>%
    mutate(estimand = factor(estimand, 
                             levels = c("mu", "mu_1", "mu_0"),
                             ordered = TRUE),
           theta = round(theta, 4),
           tau = round(tau, 4))
  
  return(results)
}





# --- Function to compute performance metrics ---
compute_metrics <- function(estimates, estimands){
  
  # 0. Compute Out-of-bound (OOB) estimates proportion and filter out
  estimates <- estimates %>%
    dplyr::mutate(oob = (is.na(est)) | (est > 1) | (est < 0))
  
  result_oob = estimates %>%
    dplyr::group_by(estimand, theta, tau) %>%
    dplyr::summarise(OOBperc = mean(oob) * 100, .groups = "drop")
  
  estimates_filtered = estimates %>% filter(!oob)
  
  # 1. Add Delta-method based PCI 
  estimates_filtered <- estimates_filtered %>%
    mutate(
      # Transform estimate to logit scale
      logit_est = qlogis(est),
      
      # Calculate the Jacobian (derivative of logit): 1 / (p * (1-p))
      jacobian = 1 / (est * (1 - est)),
      
      # Transform SE to logit scale (Delta Method)
      # se_logit = se_orig * jacobian
      se_logit = se * jacobian,
      
      # Calculate Bounds on Logit Scale
      lower_logit = logit_est - 1.96 * se_logit,
      upper_logit = logit_est + 1.96 * se_logit,
      
      # Back-transform to probability scale [0,1]
      PCL_delta = plogis(lower_logit),
      PCU_delta = plogis(upper_logit)
    ) %>%
    # Remove temporary helper columns
    dplyr::select(-logit_est, -jacobian, -se_logit, -lower_logit, -upper_logit)
  
  # 2. Join estimands
  estimates_filtered = estimates_filtered %>%
    left_join(estimands, by = c("estimand", "theta", "tau")) %>%
    mutate(
      bias = est - truth,
      cov  = I(PCL <= truth & truth <= PCU),
      cov_delta = I(PCL_delta <= truth & truth <= PCU_delta),
      ucov = I(UCL <= truth & truth <= UCU)
    )
  
  # 3. Compute UCOV (Uniform Coverage) on valid estimates only
  result_UCov = estimates_filtered %>%
    dplyr::group_by(estimand, tau, d) %>%
    dplyr::summarize(ucov = mean(ucov, na.rm = TRUE) >= 1, .groups = "drop_last") %>%
    dplyr::summarize(UCOV = round(mean(ucov, na.rm = TRUE) * 1000) / 10, .groups = "drop")
  
  # 4. Compute other metrics on valid data
  result = estimates_filtered %>%
    dplyr::group_by(estimand, theta, tau) %>%
    dplyr::summarise(Truth = mean(truth, na.rm = TRUE),
                     Bias = mean(bias, na.rm = TRUE),
                     AbsBias = mean(abs(bias), na.rm = TRUE),
                     RMSE = sqrt(mean(bias^2, na.rm = TRUE)),
                     ASE = mean(se, na.rm = TRUE),
                     ESE = sd(est, na.rm = TRUE),
                     COV = round(mean(cov, na.rm = TRUE) * 1000) / 10,
                     COV_delta = round(mean(cov_delta, na.rm = TRUE) * 1000) / 10,
                     .groups = "drop") %>%
    left_join(result_oob, by = c("estimand", "theta", "tau")) %>%
    left_join(result_UCov, by = c("estimand", "tau"))
  
  return(result)
}

# --- Main Results Processing Loop ---
data.result = data.frame()
for(N in N_values){
  print(glue("Processing N: {N}"))
  estimands_raw <- load_results(glue("estimand/N{N}")) 
  if(is.null(estimands_raw)) next
  estimands = estimands_raw %>% 
    dplyr::group_by(theta, tau, estimand) %>%
    dplyr::summarise(truth = mean(truth), .groups = "drop")
  
  estimates <- load_results(glue("estimate/TypeB_m{m}_r100_N{N}/Rdata"))
  if(is.null(estimates)) next
  
  result <- compute_metrics(estimates, estimands) %>% mutate(N = N, m = m)
  data.result = rbind(data.result, result)
}

# --- Visualization ---
data.plot_N = data.result %>%
  tidyr::pivot_longer(cols = c("Bias", "AbsBias", "ASE", "ESE", "COV", "COV_delta", "UCOV", "OOBperc"),
                      names_to = "metric", values_to = "value") %>%
  dplyr::filter(estimand %in% c("mu", "mu_1", "mu_0")) %>%
  dplyr::mutate(
    N = as.factor(N),
    estimand = factor(estimand, levels = c("mu", "mu_1", "mu_0"))
  ) %>%
  dplyr::filter(theta %in% c(0.3, 0.4, 0.5, 0.6, 0.7), tau %in% c(0.2, 0.4)) %>%
  dplyr::mutate(theta = paste0("alpha == ", theta), tau = glue("tau == {tau}"))

ggplot(data.plot_N %>% 
         dplyr::filter(metric %in% c("OOBperc", "Bias", "ESE", "COV_delta")) %>% 
         dplyr::mutate(metric = ifelse(metric == "COV_delta", "COV", metric)) %>%
         dplyr::mutate(metric = factor(metric, levels = c("OOBperc", "Bias", "ESE", "COV")))
       , 
       aes(x = estimand, y = value, fill = N)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  
  ### 95% CI benchmark line ###
  geom_hline(data = expand.grid(metric = "COV"), 
             aes(yintercept = 95),
             color = "black",
             linetype = "dashed") +
  
  scale_x_discrete(labels = c("mu" = expression(mu), "mu_1" = expression(mu[1]), "mu_0" = expression(mu[0]))) + 
  facet_nested(metric ~ theta, scales = "free", labeller = label_parsed) +
  
  # Force OOBperc to 0-100 scale
  facetted_pos_scales(
    y = list(metric == "OOBperc" ~ scale_y_continuous(limits = c(0, 100)))
  ) +
  
  theme_bw() + 
  theme(axis.title.y = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "gray95")) +
  labs(title = "Simulation Performance Metrics (Valid Estimates Only)")

ggsave(filename = glue("FigS.Small_m{m}_Large_N{paste(N_values, collapse = ',')}.pdf"), width = 9, height = 5)




#-----------------------------------------------------------------------------#
#--- Analysis of histogram ---#
#-----------------------------------------------------------------------------#

# 1. Define the list of N values you want to analyze


# Initialize a list to store plots
plot_list <- list()

# 2. Loop through each N
for (N in N_values) {
  
  # --- Load Data ---
  # Wrap in tryCatch to skip missing files gracefully
  try({
    estimates.N <- load_results(glue("estimate/TypeB_m{m}_r100_N{N}/Rdata"))
    
    if (is.null(estimates.N)) next
    
    # --- Filter Data ---
    check_data <- estimates.N %>% 
      filter(estimand == "mu", theta == 0.5, tau == 0.2) %>%
      filter(!is.na(est)) %>%
      filter(abs(est) < 1) # Exclude extreme estimates
    
    # Skip if not enough data
    if (nrow(check_data) < 10) next
    
    # --- Calculate Statistics ---
    x <- check_data$est
    
    # 1. Shapiro-Wilk (limit to 5000 samples max)
    sw_test <- shapiro.test(x[1:min(5000, length(x))])
    sw_p <- sw_test$p.value
    
    # 2. Skewness (Sample Skewness)
    n_samples <- length(x)
    skew_val <- (sum((x - mean(x))^3) / n_samples) / (sd(x)^3)
    
    # 3. Kurtosis (Sample Kurtosis - Normal is approx 3)
    kurt_val <- (sum((x - mean(x))^4) / n_samples) / (sd(x)^4)
    
    # Create Label String
    stats_label <- glue("SW p-value: {signif(sw_p, 3)} | Skewness: {round(skew_val, 3)} | Kurtosis: {round(kurt_val, 3)}")
    
    # --- Plot 1: Density with Normal Overlay ---
    p1 <- ggplot(check_data, aes(x = est)) +
      geom_histogram(aes(y = ..density..), bins = 60, fill = "lightblue", color = "black", alpha = 0.7) +
      geom_density(alpha = 0.4, fill = "blue", size = 0.8) +
      stat_function(fun = dnorm, 
                    args = list(mean = mean(x), sd = sd(x)), 
                    color = "red", linetype = "dashed", size = 1) +
      theme_bw() +
      theme(plot.subtitle = element_text(size = 9, face = "italic")) +
      labs(title = glue("Distribution (N={N}, m={m})"),
           subtitle = stats_label,
           x = "Estimated mu", y = "Density")
    
    # --- Plot 2: Q-Q Plot ---
    p2 <- ggplot(check_data, aes(sample = est)) +
      stat_qq(size = 1, alpha = 0.5) +
      stat_qq_line(color = "red", size = 1) +
      theme_bw() +
      labs(title = glue("Q-Q Plot (N={N}, m={m})"),
           subtitle = "Reference line: Normal Distribution",
           x = "Theoretical", y = "Sample")
    
    # Add to list (P1 then P2)
    plot_list[[length(plot_list) + 1]] <- p1
    plot_list[[length(plot_list) + 1]] <- p2
    
  }, silent = TRUE)
}

# 3. Arrange and Display all plots
# ncol=2 ensures Density is on Left, QQ is on Right for each N
if (length(plot_list) > 0) {
  grid.arrange(grobs = plot_list, ncol = 2)
} else {
  print("No data found or processed.")
}