root.dir = "~/research/NPSACI"
setwd(paste0(root.dir,"/simulation/M.2.IPCW_comparison"))

###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(dplyr)

thetas <- 0.1*1:9
taus <- c(0.3)

###------ Read estimand ------###
file.list <- list.files(paste0(root.dir, "/simulation/M.1.TypeB_Acorr/estimand"), pattern = "estimand.*rds", full.names = T)

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

###------ SBS-NSS estimation result ------###

estimate.list = list.files(paste0(root.dir, "/simulation/M.1.TypeB_Acorr/estimate/m200_r100_nonpara/Rdata"), 
                           pattern = "estimate.*rds", full.names = T)

result.list <- lapply(estimate.list, function(file) readRDS(file)$result)

print(paste0(length(result.list), " estimate Rdata files were loaded"))

## Aggregate estimation results
result_NSS = bind_rows(result.list, .id = "d")

## Merge with true estimands values
result_NSS2 = result_NSS %>%
  left_join(estimands, by = c("estimand", "theta", "tau")) %>%
  mutate(bias = est - truth,
         cov = I(abs(est - truth) < 1.96*se))

## Summarize performance metrics
result_table_NSS = result_NSS2 %>%
  dplyr::filter(estimand %in% c("mu_1", "mu_0")) %>%
  dplyr::group_by(estimand, theta, tau) %>%
  dplyr::summarise(Truth = mean(truth, na.rm = T),
            Bias = mean(bias, na.rm = T),
            RMSE = sqrt(mean(bias^2, na.rm = T)),
            ASE = mean(se, na.rm = T),
            ESE = sd(est, na.rm = T),
            Cov = round(mean(cov, na.rm = T),2)*100) %>%
  dplyr::arrange(theta,
          factor(estimand, levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")), 
          tau) %>%
  dplyr::filter(tau == 0.3)
  
result_table_NSS


###------ SBS-PSS estimation result ------###

estimate.list = list.files(paste0(root.dir, "/simulation/M.1.TypeB_Acorr/estimate/m200_r100_para/Rdata"), 
                           pattern = "estimate.*rds", full.names = T)

result.list <- lapply(estimate.list, function(file) readRDS(file)$result)

print(paste0(length(result.list), " estimate Rdata files were loaded"))

## Aggregate estimation results
result_PSS = bind_rows(result.list, .id = "d")

## Merge with true estimands values
result_PSS2 = result_PSS %>%
  left_join(estimands, by = c("estimand", "theta", "tau")) %>%
  mutate(bias = est - truth,
         cov = I(abs(est - truth) < 1.96*se))

## Summarize performance metrics
result_table_PSS = result_PSS2 %>%
  dplyr::filter(estimand %in% c("mu_1", "mu_0")) %>%
  dplyr::group_by(estimand, theta, tau) %>%
  dplyr::summarise(Truth = mean(truth, na.rm = T),
                   Bias = mean(bias, na.rm = T),
                   RMSE = sqrt(mean(bias^2, na.rm = T)),
                   ASE = mean(se, na.rm = T),
                   ESE = sd(est, na.rm = T),
                   Cov = round(mean(cov, na.rm = T),2)*100) %>%
  dplyr::arrange(theta,
                 factor(estimand, levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")), 
                 tau) %>%
  dplyr::filter(tau == 0.3)

result_table_PSS


###------ Chakladar IPCW estimation result ------###

estimate.list = list.files(paste0("estimate_chak/m200/Rdata"), 
                           pattern = "estimate.*rds", full.names = T)

result.list <- lapply(estimate.list, function(file) readRDS(file))

print(paste0(length(result.list), " estimate Rdata files were loaded"))

## Change result format to long format

for(d in 1:length(result.list)){
  result_sub = as.data.frame(result.list[[d]]) %>% 
    tibble::rownames_to_column("theta") %>%
    dplyr::mutate(theta = as.numeric(theta))
  
  result_sub1 = result_sub %>% 
    dplyr::select(theta, "mu-hat1", ASE1) %>% 
    dplyr::rename(est = "mu-hat1", se = ASE1) %>% 
    dplyr::mutate(estimand = "mu_1")
  
  result_sub0 = result_sub %>% 
    dplyr::select(theta, "mu-hat0", ASE0) %>% 
    dplyr::rename(est = "mu-hat0", se = ASE0) %>% 
    dplyr::mutate(estimand = "mu_0")
  
  result_sub10 = rbind(result_sub1, result_sub0)
  
  result.list[[d]] = result_sub10
}

## Aggregate estimation results
result_chak = bind_rows(result.list, .id = "d")

## Merge with true estimands values
result_chak2 = result_chak %>%
  left_join(estimands, by = c("estimand", "theta")) %>%
  mutate(bias = est - truth,
         cov = I(abs(est - truth) < 1.96*se))

## Summarize performance metrics
result_table_chak = result_chak2 %>%
  dplyr::group_by(estimand, theta, tau) %>%
  dplyr::summarise(Truth = mean(truth),
            Bias = mean(bias),
            RMSE = sqrt(mean(bias^2)),
            ASE = mean(se),
            ESE = sd(est),
            Cov = round(mean(cov, na.rm = T),2)*100) %>%
  dplyr::arrange(theta,
          factor(estimand, levels = c("mu", "mu_1", "mu_0", "de", "se_1", "se_0", "oe", "te")), 
          tau) %>%
  dplyr::filter(tau == 0.3)

result_table_chak



### Comparison ###
result_table_combined = inner_join(result_table_NSS, 
                                   result_table_PSS, 
                                   by = c("estimand", "theta", "tau", "Truth"), 
                                   suffix = c(".NSS", ".PSS"))

result_table_combined = inner_join(result_table_combined, 
                                   result_table_chak, 
                                   by = c("estimand", "theta", "tau", "Truth"), 
                                   suffix = c("", ".chak"))

# Identify continuous columns (numeric or integer)
continuous_cols <- sapply(result_table_combined, is.numeric)

# Round up to 3 decimals for continuous columns
result_table_combined[continuous_cols] <- round(result_table_combined[continuous_cols], 3)

result_table_combined


### Visulaization ###
data.plot = rbind(
  result_table_NSS %>% dplyr::mutate(method = "NSS"),
  result_table_PSS %>% dplyr::mutate(method = "PSS"),
  result_table_chak %>% dplyr::mutate(method = "Chakladar")
) %>% select(method, theta, estimand, Bias, RMSE, Cov) %>% ungroup()


library(reshape)
data.plot = as.data.frame(data.plot)
data.plot.2 = reshape::melt(data.plot, id.vars = c("method", "theta", "estimand"), var = "type") %>%
  filter(theta %in% c(0.1,0.3,0.5,0.7,0.9)) %>%
  mutate(method = factor(method, levels = c("NSS", "PSS", "Chakladar"), ordered = TRUE))

### Bar chart (NSS vs P IPW) ###
ggplot(data = data.plot.2 %>% 
         mutate(theta = paste0("alpha == ", theta)),
       aes(x = estimand, y = value, group = method, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ### 95% CI benchmark line ###
  geom_hline(data = expand.grid(method = c("NSS", "PSS", "Chakladar"),
                                type = "Cov",
                                theta = paste0("alpha == ", c(0.1,0.3,0.5,0.7,0.9))), 
             aes(yintercept = 95),
             color = "blue",
             linetype = "dashed") +
  facet_grid(type ~ theta, scales = "free_y", labeller = label_parsed) +
  scale_fill_discrete("Method", labels = c("SBS-NSS", "SBS-PSS", "Chakladar IPCW")) + 
  scale_x_discrete(labels = c("mu_1" = expression(mu[1](alpha)),
                              "mu_0" = expression(mu[0](alpha)))) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

ggsave(filename = "Fig1.NSS.PSS.Chakladar.comparison.barchart.pdf", width = 10, height = 5)
