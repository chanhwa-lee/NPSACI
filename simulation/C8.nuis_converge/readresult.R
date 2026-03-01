###------------------- Load libraries ----------------------###
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.5")
library(dplyr)
library(ggplot2)
library(glue)
library(conflicted)
library(tidyr)

# Set conflict preferences
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::mutate)

#--- Nuisance estimators convergence plot ---#
# m_seq = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
m_seq = c(50, 100, 200, 400, 800, 1600)

nuis_est.df = data.frame()

for(m in m_seq){
  nuis_est.file.list = list.files(glue("nuis_estimates/m{m}/Rdata"),
                                  pattern = "nuis_estimate.*rds", full.names = T)

  if(length(nuis_est.file.list) < 10){next}

  nuis_est.list <- lapply(nuis_est.file.list, function(file) readRDS(file))

  print(glue("m = {m}, {length(nuis_est.file.list)} files were loaded"))

  nuis_est = bind_rows(nuis_est.list, .id = "d")

  nuis_est.df = rbind(nuis_est.df, nuis_est %>% mutate(m = m))
}

ggplot(nuis_est.df, aes(x = factor(m), y = bias, col = method)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  theme_bw() +
  labs(
    x = "Number of clusters (m)",
    y = expression(L[2](P) ~ bias),
    col = "Method"
  ) +
  facet_wrap(
    ~nuis,
    labeller = as_labeller(
      c("F" = "bold(F)^T ~ '(Event time)'",
        "G" = "bold(S)^C ~ '(Censoring time)'",
        "H" = "H ~ '(Treatment)'"), 
      default = label_parsed),
    scales = "free_y") +
  ylim(0, NA) +
  theme(legend.position = c(0.01, 0.30),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "black"))

ggsave(filename = glue("FigS.NuisanceEstimatorsBias.pdf"), width = 12, height = 4)






###--------------- Supp: Convergence rate fitting line ------------------###

# #--- 1. Simple Linear Regression on Log-Log Scale ---
# # We calculate alpha (slope) and C (exp of intercept) for each group
# fit_stats <- nuis_est.df %>%
#   group_by(nuis, method) %>%
#   summarize(
#     model = list(lm(log(bias) ~ log(m))),
#     alpha = coef(model[[1]])[2],
#     C = exp(coef(model[[1]])[1]),
#     .groups = "drop"
#   )
# 
# # Create the points for the fitting lines using the estimated parameters
# fit_lines <- fit_stats %>%
#   rowwise() %>%
#   reframe(data.frame(
#     nuis = nuis,
#     method = method,
#     m = m_seq,
#     fit_y = C * (m_seq^alpha)
#   ))
# 
# # Create the mathematical labels for alpha
# alpha_labels <- fit_stats %>%
#   mutate(label = paste0("alpha[", method, "] == ", round(alpha, 3)))
# 
# #--- 2. Generate the Final Plot ---
# ggplot(nuis_est.df, aes(x = factor(m), y = bias, col = method)) +
#   geom_boxplot(position = position_dodge(width = 0.75)) +
#   # Add the fitting lines using the calculated fit_y
#   geom_line(data = fit_lines, aes(x = factor(m), y = fit_y, group = method, col = method), 
#             linetype = "dashed") +
#   # Add alpha labels stacked vertically
#   geom_text(data = alpha_labels, 
#             aes(x = Inf, y = Inf, label = label, col = method), 
#             hjust = 1.1, 
#             vjust = ifelse(alpha_labels$method == "para", 1.5, 3.2), 
#             parse = TRUE, show.legend = FALSE) +
#   theme_bw() +
#   labs(
#     x = "Number of clusters (m)",
#     y = expression(L[2](P) ~ bias),
#     col = "Method"
#   ) +
#   facet_wrap(
#     ~nuis,
#     labeller = as_labeller(
#       c("F" = "bold(F)^T ~ '(Event time)'",
#         "G" = "bold(S)^C ~ '(Censoring time)'",
#         "H" = "H ~ '(Treatment)'"), 
#       default = label_parsed),
#     scales = "free_y") +
#   ylim(0, NA) +
#   theme(
#     legend.position = "inside", 
#     legend.position.inside = c(0.05, 0.25),
#     legend.justification = c(0, 1),
#     legend.background = element_rect(fill = "white", color = "black")
#   )
# 
# ggsave(filename = glue("FigS.NuisanceEstimatorsBias_withFittedLine.pdf"), width = 12, height = 4)