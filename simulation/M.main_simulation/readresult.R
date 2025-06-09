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


# --- Load estimands and aggregate ---
estimands_TypeB <- load_results("estimand/TypeB") %>% 
  dplyr::group_by(theta, tau, estimand) %>%
  dplyr::summarise(truth = mean(truth), .groups = "drop")

estimands_TPB   <- load_results("estimand/TPB") %>%
  dplyr::group_by(theta, tau, estimand) %>%
  dplyr::summarise(truth = mean(truth), .groups = "drop")


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
    dplyr::summarize(ucov = mean(ucov, na.rm = TRUE) >= 0.95, .groups = "drop_last") %>%
    dplyr::summarize(UCOV = ceiling(mean(ucov, na.rm = T)*100), .groups = "drop")
  
  result = estimates %>%
    dplyr::group_by(estimand, theta, tau) %>%
    dplyr::summarise(Truth = mean(truth, na.rm = T),
                     Bias = mean(bias, na.rm = T),
                     AbsBias = mean(abs(bias), na.rm = T),
                     RelBias = abs(mean(bias, na.rm = T) / Truth),
                     RMSE = sqrt(mean(bias^2, na.rm = T)),
                     ASE = mean(se, na.rm = T),
                     ESE = sd(est, na.rm = T),
                     COV = ceiling(mean(cov, na.rm = T)*100),
                     .groups = "drop") %>%
    left_join(result_UCov, by = c("estimand", "tau"))
  
  return(result)
  
}

# --- Load estimates and combine with estimands to compute metrics ---
m = 200
r = 100
Scenarios = c(1,21,22,23,3,41,42,43)

result_TypeB = lapply(Scenarios,
                         function(s){
                           estimates = load_results(glue("estimate/TypeB_m{m}_r{r}_Scenario{s}/Rdata"))
                           if(is.null(estimates)) return()
                           compute_metrics(estimates, estimands_TypeB)  %>% mutate(Scenario = !!s)
                          }) %>% bind_rows()

result_TPB = lapply(Scenarios,
                      function(s){
                        estimates = load_results(glue("estimate/TPB_m{m}_r{r}_Scenario{s}/Rdata"))
                        if(is.null(estimates)) return()
                        compute_metrics(estimates, estimands_TPB)  %>% mutate(Scenario = !!s)
                      }) %>% bind_rows()


###------------------------------------------------------###
###------------- Table3. TypeB NCF result ---------------###
###------------------------------------------------------###

table_TypeB = result_TypeB %>% 
  filter(Scenario == 1) %>%
  filter(tau %in% c(0.2,0.4)) %>%
  filter( (theta == 0.3) |
            (theta == 0.45 & estimand %in% c("mu", "mu_1", "mu_0", "de")) |
            (theta == 0.6) ) %>%
  mutate(across(c("Truth", "Bias", "RMSE", "ASE", "ESE"), 
                ~sprintf("%.1f", ifelse(abs(.x * 100) < 1e-8, 0, .x * 100)))) %>%
  select(estimand, theta, tau, Truth, Bias, RMSE, ASE, ESE, COV, UCOV, Scenario) %>%
  arrange(tau, theta) %>%
  print(n=100)

write.csv(table_TypeB, "Table3.Simulation_TypeB.csv")

#--- Supplementary: TypeB PCF result ---
result_TypeB %>% 
  filter(Scenario == 43) %>%
  filter(tau %in% c(0.2,0.4)) %>%
  filter( (theta == 0.3) |
            (theta == 0.45 & estimand %in% c("mu", "mu_1", "mu_0", "de")) |
            (theta == 0.6) )%>%
  arrange(tau, theta) %>%
  print(n=100)


###------------------------------------------------------###
###------------- TableS1. TPPB NCF result ---------------###
###------------------------------------------------------###

table_TPB = result_TPB %>%
  filter(Scenario == 1) %>%
  filter(tau %in% c(0.2,0.4)) %>%
  filter( (theta == 0 & estimand %in% c("mu", "mu_1", "mu_0", "de")) |
            (theta == 0.25) |
            (theta == 0.5) )%>%
  mutate(across(c("Truth", "Bias", "RMSE", "ASE", "ESE"), 
                ~sprintf("%.1f", ifelse(abs(.x * 100) < 1e-8, 0, .x * 100)))) %>%
  select(estimand, theta, tau, Truth, Bias, RMSE, ASE, ESE, COV, UCOV, Scenario) %>%
  arrange(tau, theta) %>%
  print(n=100)

write.csv(table_TPB, "TableS1.Simulation_TPB.csv")




###-------------------------------------------------------###
###---- Fig1. NCF vs PCF vs Chakladar IPCW comparison ----###
###-------------------------------------------------------###

estimates_TypeB_chak = rbind(
  load_results(glue("estimate/TypeB_m200_Chakladar_tau0.2/Rdata")),
  load_results(glue("estimate/TypeB_m200_Chakladar_tau0.4/Rdata")))

result_TypeB_chak = compute_metrics(estimates_TypeB_chak, estimands_TypeB) %>% mutate(Scenario = "Chakladar")

result_TypeB_method = rbind(result_TypeB %>% 
                                  filter(Scenario %in% c(1,43),
                                         estimand %in% c("mu_1", "mu_0"),
                                         theta %in% c(0.3,0.45,0.6),
                                         tau %in% c(0.2,0.4)),
                                result_TypeB_chak)

# Long form to generate plotting dataset
data.plot_method =
  reshape::melt(as.data.frame(result_TypeB_method),
                id.vars = c("estimand", "theta", "tau", "Scenario"),
                var = "metric") %>%
  dplyr::filter(metric %in% c("Bias", "ESE", "COV"))

label_estimands = c("mu"   = expression(mu),
                    "mu_1" = expression(mu[1]),
                    "mu_0" = expression(mu[0]))

ggplot(data.plot_method %>%
         dplyr::mutate(theta = paste0( "alpha == ", theta),
                       tau   = glue("tau == {tau}")),
       aes(x = estimand, y = value, group = Scenario, fill = Scenario))+

  geom_bar(stat = "identity", position = position_dodge()) +

  ### 95% CI benchmark line ###
  geom_hline(data = expand.grid(metric = "COV"),
             aes(yintercept = 95),
             color = "black",
             linetype = "dashed") +

  scale_fill_discrete("Method",
                      labels = c("1"  = "SBS-NCF",
                                 "43" = "SBS-PCF",
                                 "Chakladar" = "Chakladar IPCW")) +

  scale_x_discrete(labels = label_estimands) +
  labs(x = "Estimand") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  ggh4x::facet_nested(metric ~ theta+ tau, scales = "free", labeller = label_parsed)

ggsave(filename = glue("Fig1.NCF.PCF.IPCWcomparison.TypeB.pdf"), width = 8, height = 4)



#--- Supplementary: Multiple robustness ---

result_TypeB_method_all = rbind(result_TypeB %>% 
                              filter(Scenario %in% c(1,23,3,43),
                                     estimand %in% c("mu_1", "mu_0"),
                                     theta %in% c(0.3,0.45,0.6),
                                     tau %in% c(0.2,0.4)),
                            result_TypeB_chak)

# Long form to generate plotting dataset
data.plot_method_all =
  reshape::melt(as.data.frame(result_TypeB_method_all),
                id.vars = c("estimand", "theta", "tau", "Scenario"),
                var = "metric") %>%
  dplyr::filter(metric %in% c("Bias", "ESE", "COV"))

ggplot(data.plot_method_all %>%
         dplyr::mutate(theta = paste0( "alpha == ", theta),
                       tau   = glue("tau == {tau}")),
       aes(x = estimand, y = value, group = Scenario, fill = Scenario))+
  
  geom_bar(stat = "identity", position = position_dodge()) +
  
  ### 95% CI benchmark line ###
  geom_hline(data = expand.grid(metric = "COV"),
             aes(yintercept = 95),
             color = "black",
             linetype = "dashed") +
  
  scale_fill_discrete("Method",
                      labels = c("1"  = "SBS-NCF",
                                 "43" = "SBS-PCF",
                                 "Chakladar" = "Chakladar IPCW")) +
  
  scale_fill_discrete("Method",
                      labels = c("1"  = "CF: T np / C np / A np",
                                 "23" = "CF: T np / C p  / A p",
                                 "3"  = "CF: T p  / C np / A np",
                                 "43" = "CF: T p  / C p  / A p",
                                 "Chakladar" = "Chakladar IPCW")) +
  
  scale_x_discrete(labels = label_estimands) +
  labs(x = "Estimand") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  ggh4x::facet_nested(metric ~ theta+ tau, scales = "free", labeller = label_parsed)




###-------------------------------------------------------###
###--------------- FigS1. Performance by m ---------------###
###-------------------------------------------------------###

result_TypeB_m = lapply(c(25,50,100,200,400),
                        function(m){
                          estimates = load_results(glue("estimate/TypeB_m{m}_r100_Scenario1/Rdata"))
                          if(is.null(estimates)) return()
                          compute_metrics(estimates, estimands_TypeB) %>% mutate(m = !!m)
                        }) %>% bind_rows()

# Long form to generate plotting dataset
data.plot_m =
  reshape::melt(as.data.frame(result_TypeB_m),
                id.vars = c("estimand", "theta", "tau", "m"),
                var = "metric") %>%
  dplyr::filter(metric %in% c("Bias", "ESE", "COV", "UCOV"),
                estimand %in% c("mu", "mu_1", "mu_0")) %>%
  dplyr::mutate(m = factor(m))

ggplot(data.plot_m %>%
         dplyr::filter(theta %in% c(0.3,0.45,0.6),
                       tau   %in% c(0.2,0.4)) %>%
         dplyr::mutate(theta = paste0( "alpha == ", theta),
                       tau   = glue("tau == {tau}")),
       aes(x = estimand, y = value, group = m, fill = m))+
  
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
  
  # scale_fill_discrete("m") +
  scale_x_discrete(labels = label_estimands) +
  labs(x = "Estimand") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  ggh4x::facet_nested(metric ~ theta+ tau, scales = "free", labeller = label_parsed)

ggsave(filename = glue("FigS1.SimulationOver_m.TypeB.pdf"), width = 9, height = 5)



###-------------------------------------------------------###
###------------ FigS2. Performance by bounding -----------###
###-------------------------------------------------------###

estimates_TypeB_unbounded = load_results(glue("estimate/TypeB_m200_r100_Scenario1_unbounded/Rdata"))
result_TypeB_unbounded = compute_metrics(estimates_TypeB_unbounded, estimands_TypeB) %>% mutate(Scenario = "Unbounded")
  
result_TypeB_bounding = rbind(result_TypeB %>% filter(Scenario == 1) %>% mutate(Scenario = "Bounded"),
                              result_TypeB_unbounded)

result_TypeB_bounding

# Long form to generate plotting dataset
data.plot_bounding =
  reshape::melt(as.data.frame(result_TypeB_bounding),
                id.vars = c("estimand", "theta", "tau", "Scenario"),
                var = "metric") %>%
  dplyr::filter(metric %in% c("Bias", "ESE", "COV", "UCOV"),
                estimand %in% c("mu", "mu_1", "mu_0"))

ggplot(data.plot_bounding %>%
         dplyr::filter(theta %in% c(0.3,0.45,0.6),
                       tau   %in% c(0.2,0.4)) %>%
         dplyr::mutate(theta = paste0( "alpha == ", theta),
                       tau   = glue("tau == {tau}")),
       aes(x = estimand, y = value, group = Scenario, fill = Scenario))+

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

  scale_x_discrete(labels = label_estimands) +
  scale_fill_discrete(name = "Method") + 
  labs(x = "Estimand") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  ggh4x::facet_nested(metric ~ theta+ tau, scales = "free", labeller = label_parsed)

ggsave(filename = glue("FigS2.BoundingComparison.pdf"), width = 9, height = 5)




###-------------------------------------------------------###
###--------------- FigS3. Performance by r ---------------###
###-------------------------------------------------------###

result_TypeB_r = lapply(c(10,20,50,100,200,500),
                        function(r){
                          estimates = load_results(glue("estimate/TypeB_m200_r{r}_Scenario1/Rdata"))
                          if(is.null(estimates)) return()
                          compute_metrics(estimates, estimands_TypeB) %>% mutate(r = !!r)
                        }) %>% bind_rows()

# Long form to generate plotting dataset
data.plot_r =
  reshape::melt(as.data.frame(result_TypeB_r),
                id.vars = c("estimand", "theta", "tau", "r"),
                var = "metric") %>%
  dplyr::filter(metric %in% c("Bias", "ESE", "COV", "UCOV"),
                estimand %in% c("mu", "mu_1", "mu_0")) %>%
  dplyr::mutate(r = factor(r))

ggplot(data.plot_r %>%
         dplyr::filter(theta %in% c(0.3,0.45,0.6),
                       tau   %in% c(0.2,0.4)) %>%
         dplyr::mutate(theta = paste0( "alpha == ", theta),
                       tau   = glue("tau == {tau}")),
       aes(x = estimand, y = value, group = r, fill = r))+

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

  scale_fill_discrete(name = "r", guide = guide_legend(nrow = 1)) +
  scale_x_discrete(labels = label_estimands) +
  labs(x = "Estimand") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  ggh4x::facet_nested(metric ~ theta+ tau, scales = "free", labeller = label_parsed)


ggsave(filename = glue("FigS3.SimulationOver_r.TypeB.pdf"), width = 9, height = 5)









###-------------------------------------------------------###
###------ Supplementary: Performance over nuis est -------###
###-------------------------------------------------------###

make_performance_plot = function(result, thetas, taus, policy){
  
  # Long form to generate plotting dataset
  data.plot = 
    reshape::melt(as.data.frame(result), 
                  id.vars = c("estimand", "theta", "tau", "Scenario"), 
                  var = "metric") %>%
    dplyr::filter(metric %in% c("RelBias", "ASE", "ESE", "COV", "UCOV"), 
                  estimand %in% c("mu", "mu_1", "mu_0")) %>%
    dplyr::mutate(Scenario = factor(Scenario, 
                                    levels = c("1", "21", "22", "23", "3", "41", "42", "43")))
  
  ggplot(data.plot %>% 
           dplyr::filter(theta %in% thetas,
                         tau   %in% taus) %>%
           dplyr::mutate(theta = paste0( ifelse(policy == "TypeB", "alpha == ", "rho == "), theta),
                         tau   = glue("tau == {tau}")), 
         aes(x = estimand, y = value, group = Scenario, fill = as.factor(Scenario)))+ 
    
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
    
    scale_fill_discrete("Nuisance estimators",
                        labels = c("1"  = "1 : T np / C np / A np",
                                   "21" = "21: T np / C np / A p",
                                   "22" = "22: T np / C p  / A np",
                                   "23" = "23: T np / C p  / A p",
                                   "3"  = "3 : T p  / C np / A np",
                                   "41" = "41: T p  / C np / A p",
                                   "42" = "42: T p  / C p  / A np",
                                   "43" = "43: T p  / C p  / A p")) +
    scale_x_discrete(labels = label_estimands) + 
    labs(x = "Estimand") + 
    theme_bw() + 
    theme(axis.title.y = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal") + 
    ggh4x::facet_nested(metric ~ theta+ tau, scales = "free", labeller = label_parsed)
}


make_performance_plot(result_TypeB, thetas = c(0.3,0.45,0.6), taus = c(0.2,0.4), "TypeB")
ggsave(filename = glue("FigS.SimulationOverNuisanceEstimators.TypeB.m{m}.r{r}.pdf"), width = 12, height = 8)

make_performance_plot(result_TPB, thetas = c(0, 0.025*10, 0.5), taus = c(0.2,0.4), "TPB")
ggsave(filename = glue("FigS.SimulationOverNuisanceEstimators.TPB.m{m}.r{r}.pdf"), width = 12, height = 8)


make_performance_plot(result_TypeB %>% filter(Scenario %in% c(1,23,3,43)), thetas = c(0.3,0.45,0.6), taus = c(0.2,0.4), "TypeB")






# #--- Supplementary: Nuisance estimators convergence plot ---
# nuis_est.df = data.frame()
# 
# for(m in c(25, 50, 100, 200, 400, 800)){
#   nuis_est.file.list = list.files(glue("nuis_estimates/m{m}/Rdata"),
#                                   pattern = "nuis_estimate.*rds", full.names = T)
# 
#   if(length(nuis_est.file.list) < 10){next}
# 
#   nuis_est.list <- lapply(nuis_est.file.list, function(file) readRDS(file))
# 
#   print(glue("m = {m}, {length(nuis_est.file.list)} files were loaded"))
# 
#   nuis_est = bind_rows(nuis_est.list, .id = "d")
# 
#   nuis_est.df = rbind(nuis_est.df, nuis_est %>% mutate(m = m))
# }
# 
# ggplot(nuis_est.df, aes(x = factor(m), y = bias, col = method)) +
#   geom_boxplot(position = position_dodge(width = 0.75)) +
#   theme_bw() +
#   labs(
#     x = "Number of clusters (m)",
#     y = "Average Relative Bias",
#     col = "Method"
#   ) +
#   facet_wrap(~nuis) +
#   ylim(c(0,1)) +
#   theme(legend.position = c(0.05, 0.95),  # (x, y) in [0,1] coordinates: left (0.05), top (0.95)
#         legend.justification = c(0,1),    # anchor the legend box to top-left
#         legend.background = element_rect(fill = "white", color = "black"))  # white background box
# 
# ggsave(filename = glue("FigS.NuisanceEstimatorsBias.pdf"), width = 12, height = 4)



# --- Supplementary: Plot Pointwise CI and Uniform confidence bound for one estimation result ---

# policy = "TypeB"
# m = 200
# r = 100
# s = 1 ## Nuisance estimation scenario
# d = 1 ## Simulation idx
# 
# drawCIandUCB = function(policy, estimands, m, r, s, d){
#   
#   result <- readRDS(glue("estimate/{policy}_m{m}_r{r}_Scenario{s}/Rdata/estimate_id{d}.rds")) %>%
#     ## Ad-hoc solution for floating point issue in theta and tau values (0.1*3 != 0.3)
#     dplyr::mutate(tau = round(tau,4),
#                   theta = round(theta,4))
#   
#   if (policy == "TypeB") {
#     label_estimands = c("mu"   = expression(mu[B](tau:alpha)),
#                         "mu_1" = expression(mu[B*","*1](tau:alpha)),
#                         "mu_0" = expression(mu[B*","*0](tau:alpha)),
#                         "de" = expression(DE[B](tau:alpha)),  
#                         "se_1" = expression(SE[B*","*1](tau:alpha, alpha[0])), 
#                         "se_0" = expression(SE[B*","*0](tau:alpha, alpha[0])), 
#                         "oe" = expression(OE[B](tau:alpha, alpha[0])))
#     
#   } else if (policy == "TPB") {
#     label_estimands = c("mu"   = expression(mu[TPB](tau:rho)),
#                         "mu_1" = expression(mu[TPB*","*1](tau:rho)),
#                         "mu_0" = expression(mu[TPB*","*0](tau:rho)),
#                         "de" = expression(DE[TPB](tau:rho)),  
#                         "se_1" = expression(SE[TPB*","*1](tau:rho, rho[0])), 
#                         "se_0" = expression(SE[TPB*","*0](tau:rho, rho[0])), 
#                         "oe" = expression(OE[TPB](tau:rho, rho[0])))
#   } else {
#     stop("Unknown policy type")
#   }
#   
#   ## Merge with true estimands values
#   result = result %>%
#     left_join(estimands, by = c("estimand", "theta", "tau")) %>%
#     dplyr::mutate(estimand = factor(estimand, 
#                                     levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe"),
#                                     labels = label_estimands)) 
#   
#   scenario = c("1"  = "1 : T np / C np / A np",
#                "21" = "21: T np / C np / A p",
#                "22" = "22: T np / C p  / A np",
#                "23" = "23: T np / C p  / A p",
#                "3"  = "3 : T p  / C np / A np",
#                "41" = "41: T p  / C np / A p",
#                "42" = "42: T p  / C p  / A np",
#                "43" = "43: T p  / C p  / A p")
#   
#   ggplot(result %>%
#            dplyr::filter(tau %in% c(0.1,0.3,0.5)) %>%
#            dplyr::mutate(tau = glue("tau == {tau}")), 
#          aes(x = theta, y = est)) +
#     ## Uniform confidence band
#     geom_ribbon(aes(ymin = UCL, ymax = UCU, fill = "Uniform 95% CI"), alpha = 1) +
#     
#     ## Pointwise confidence interval
#     geom_ribbon(aes(ymin = PCL, ymax = PCU, fill = "Pointwise 95% CI"), alpha = 1) +
#     
#     ## Estimate line
#     geom_line(aes(color = "Estimate"), linewidth = 0.5) +
#     
#     ## True value line
#     geom_line(aes(y = truth, color = "Truth"), linewidth = 0.5, linetype = "dashed") +
#     
#     ## Facet
#     facet_grid(estimand ~ tau, scales = "free", labeller = label_parsed) +
#     
#     ## Manual legend control
#     scale_fill_manual(name = NULL,
#                       values = c("Pointwise 95% CI" = "lightblue3",
#                                  "Uniform 95% CI" = "lightblue1")) +
#     scale_color_manual(name = NULL,
#                        values = c("Estimate" = "black",
#                                   "Truth" = "red")) +
#     
#     ## Labels and theme
#     labs(x = ifelse(policy == "TypeB", expression(alpha), expression(rho)),
#          y = "Estimate",
#          title = "Estimates with Pointwise and Uniform Confidence Bands",
#          subtitle = glue("m = {m}, simul.id = {d}, nuis.est.scenario = {scenario[as.character(s)]}")) +
#     theme_bw() +
#     theme(legend.position = "bottom",
#           legend.background = element_rect(color = "black"),
#           legend.key = element_blank())
# }
# 
# plot_grid(drawCIandUCB(policy="TypeB",estimands_TypeB,m=200,r=100,s=1,d=1),
#           drawCIandUCB(policy="TypeB",estimands_TypeB,m=200,r=100,s=43,d=1),
#           nrow = 1)
# 
# ggsave("Fig.UCB_simulation_example_TypeB.pdf", width = 12, height = 8)

# plot_grid(drawCIandUCB(policy="TPB",estimands_TPB,m=200,r=100,s=1,d=1),
#           drawCIandUCB(policy="TPB",estimands_TPB,m=200,r=100,s=22,d=1),
#           nrow = 1)
# 
# ggsave("Fig.UCB_simulation_example_TPB.pdf", width = 12, height = 8)




