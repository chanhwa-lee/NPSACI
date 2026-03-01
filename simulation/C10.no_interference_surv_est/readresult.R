###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.5")
library(dplyr)
library(tidyr)     # Replaces reshape
library(ggplot2)
library(conflicted)
library(glue)
library(cowplot)
library(latex2exp)
library(purrr)     # For map_dfr

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::mutate)

###------------------- Helper Functions ----------------------###

# 1. Generic result loader
load_results <- function(path) {
  files <- list.files(path, pattern = ".*rds", full.names = TRUE)
  if (length(files) < 10) return(NULL)
  message(length(files), " Rdata files loaded from ", path)
  
  map_dfr(files, readRDS, .id = "d") %>%
    mutate(estimand = factor(estimand, levels = c("mu", "mu_1", "mu_0", "de", "oe", "se_1", "se_0"), ordered = TRUE),
           theta = round(theta, 4),
           tau = round(tau, 4))
}

# 2. Unified Performance Metric Calculator
calc_metrics <- function(data, truth_df) {
  data %>%
    left_join(truth_df, by = c("estimand", "theta", "tau")) %>%
    mutate(bias = est - truth,
           cov  = I(PCL <= truth & truth <= PCU)) %>%
    group_by(estimand, theta, tau, method) %>%
    summarise(Truth   = mean(truth, na.rm = TRUE),
              Bias    = mean(bias, na.rm = TRUE),
              AbsBias = mean(abs(bias), na.rm = TRUE),
              RelBias = abs(mean(bias, na.rm = TRUE) / Truth),
              RMSE    = sqrt(mean(bias^2, na.rm = TRUE)),
              ASE     = mean(se, na.rm = TRUE),
              ESE     = sd(est, na.rm = TRUE),
              COV     = round(mean(cov, na.rm = TRUE) * 1000) / 10,
              .groups = "drop")
}

# 3. Plot Data Prepper (Pivots and sets factors)
prep_plot_data <- function(result_df) {
  result_df %>%
    pivot_longer(cols = c(Bias, ESE, COV), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c("Bias", "ESE", "COV")),
           method = factor(method))
}

# 4. Helper for the dummy limits trick
get_dummy_limits <- function(plot_df) {
  data.frame(
    metric = factor(c("COV", "COV", "ESE"), levels = c("Bias", "ESE", "COV")),
    value  = c(0, 100, 0),
    method = unique(plot_df$method)[1],
    theta  = if("theta" %in% names(plot_df)) unique(plot_df$theta)[1] else 1
  )
}

# 5. Shared Theme
common_theme <- theme_bw() + 
  theme(axis.title.y = element_blank())


###------------------- A. Direct Effects (DE) ----------------------###

# --- Load Data ---
estimands_de <- load_results("estimand/TypeB_alpha_0.3_0.6") %>% 
  group_by(theta, tau, estimand) %>%
  summarise(truth = mean(truth), .groups = "drop") %>%
  filter(tau == 0.2, estimand == "de")

est_prop_de <- load_results(glue("estimate/TypeB_m200_r100_alpha_0.3_0.6/Rdata")) %>%
  filter(tau == 0.2, estimand == "de") %>%
  select(-UCL, -UCU) %>%
  mutate(method = "Proposed")

est_ni_de <- map_dfr(list.files("estimate/No_interference_m200/Rdata", full.names = T), readRDS, .id = "d") %>%
  filter(Method != "AIPTW") %>%
  mutate(estimand = "de", tau = 0.2, est = -Estimate, 
         PCL = -Upper_95_CI, PCU = -Lower_95_CI, 
         se = (PCU - PCL) / (2 * 1.96), method = Method) %>%
  cross_join(expand.grid(theta = round(0.01 * (30:60), 4))) %>%
  select(d, estimand, theta, tau, est, se, PCL, PCU, method)

# --- Process ---
res_de <- bind_rows(est_prop_de, est_ni_de) %>% 
  calc_metrics(truth_df = estimands_de)

# --- Plot ---
plot_de_data <- prep_plot_data(res_de)
dummy_de <- get_dummy_limits(plot_de_data)

p.DE <- ggplot(plot_de_data, aes(x = theta, y = value, group = method, color = method)) +
  geom_point(aes(shape = method)) +
  geom_line() +
  geom_hline(data = subset(plot_de_data, metric == "Bias"), aes(yintercept = 0), color = "red", linetype = "dashed") +
  geom_hline(data = subset(plot_de_data, metric == "COV"), aes(yintercept = 95), color = "black", linetype = "dashed") +
  geom_blank(data = dummy_de) +
  labs(x = "alpha", title = "A. Direct Effect Estimation") +
  common_theme +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.direction = "horizontal",
    legend.margin = margin(4, 4, 4, 4)
  ) +
  facet_grid(metric ~ ., scales = "free")


###------------------- B. Overall Effects (OE) ----------------------###

# --- Load Data ---
estimands_oe <- load_results("estimand/TypeB_alpha_0_1") %>% 
  group_by(theta, tau, estimand) %>%
  summarise(truth = mean(truth), .groups = "drop") %>%
  filter(tau == 0.2, estimand == "oe", theta == 1)

est_prop_oe <- load_results(glue("estimate/TypeB_m200_r100_alpha_0_1/Rdata")) %>%
  filter(tau == 0.2, estimand == "oe", theta == 1) %>%
  select(-UCL, -UCU) %>%
  mutate(method = "Proposed")

est_ni_oe <- map_dfr(list.files("estimate/No_interference_m200/Rdata", full.names = T), readRDS, .id = "d") %>%
  mutate(estimand = "oe", theta = 1, tau = 0.2,
         est = ifelse(Method == "AIPTW", Estimate, -Estimate),
         PCL = ifelse(Method == "AIPTW", Lower_95_CI, -Upper_95_CI),
         PCU = ifelse(Method == "AIPTW", Upper_95_CI, -Lower_95_CI),
         se = (PCU - PCL) / (2 * 1.96), method = Method) %>%
  select(d, estimand, theta, tau, est, se, PCL, PCU, method)

# --- Process ---
res_oe <- bind_rows(est_prop_oe, est_ni_oe) %>% 
  calc_metrics(truth_df = estimands_oe)

# --- Plot ---
plot_oe_data <- prep_plot_data(res_oe %>% filter(method != "AIPTW"))
dummy_oe <- get_dummy_limits(plot_oe_data)

p.OE <- ggplot(plot_oe_data, aes(x = method, y = value, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.3) +
  geom_hline(data = subset(plot_oe_data, metric == "Bias"), aes(yintercept = 0), color = "red", linetype = "dashed") +
  geom_hline(data = subset(plot_oe_data, metric == "COV"), aes(yintercept = 95), color = "black", linetype = "dashed") +
  geom_blank(data = dummy_oe) +
  labs(x = "Method", title = "B. Overall Effect Estimation") +
  common_theme +
  theme(legend.position = "none") +
  facet_grid(metric ~ ., scales = "free")

###------------------- Combine & Save ----------------------###
plot_grid(p.DE, p.OE, nrow = 1, rel_widths = c(0.65, 0.35))
ggsave(filename = glue("FigS.Simulation_No_Interference_Method.TypeB.pdf"), width = 11, height = 5)
