library(dplyr)
library(glue)
library(ggplot2)
library(doMC)
library(randomForestSRC)
library(dbarts)

set.seed(1234)

# --- Generate synthetic data ---
generate_cluster_data = function(N){
  
  ## Step1. Covariate generation
  age <- round(runif(N, 15, 65)); dist.river <- runif(N, 0, 5)
  
  ## Step2. Treatment model
  b = rnorm(1,0,0.5)
  pi = plogis(0.2 + 0.2*(age/40-1)^2 + 0.2*pmax(dist.river/5,0.3) + b)
  A <- rbinom(N, 1, pi)
  g.A <- (sum(A)-A)/(N-1)
  
  ## Step3. Survival time model
  T_ <- round(100*rexp(N, 1/exp(2 + 1*A + 1*g.A + 1*A*g.A + 0.5*dist.river - 0.5*(age/40-1))))
  C <- round(runif(N, 0, 500 + 200*A + 50*dist.river + 100*(age/40-1)))
  
  ## Step 4. Combine all data
  Y <- pmin(T_,C)
  D <- as.numeric(T_ <= C)
  
  return(data.frame(Y = Y, D = D, A = A, age = age, dist.river = dist.river))
  
}

m = 200

toy_data = lapply(sample(x = 3:5, size = m, replace = T), generate_cluster_data) %>%
  dplyr::bind_rows(.id = "id") %>% mutate(id = as.numeric(id))


# --- Compute estimates ---

## Help functions for estimator main functions
source("code/help_util.R")

## Help functions for policy specific functions
source("code/help_TypeB.R")
source("code/help_TPB.R")

## Help functions for nuisance functions estimation method
source("code/help_nuis_est.R")

## Compute estimates
result = estimator(data = toy_data,
                   X.T.names = c("age", "dist.river"),
                   X.C.names = c("age"),
                   X.A.names = c("age", "dist.river"),
                   policy = "TypeB",
                   taus = 20*(1:20), thetas = seq(0.3, 0.6, length.out = 31), theta0 = 0.45,
                   parallel_computing = FALSE)

## Estimation Result
result$result %>% 
  filter(estimand == "mu", tau == 360) %>% 
  mutate(dplyr::across(c(est, se, PCL, PCU, UCL, UCU), ~ round(.x, 4)))

# --- Visualize estimates ---
estimates = result$result %>% mutate(theta = round(theta, 4))

## Risk over time plot
thetas = c(0.3, 0.45, 0.6)

ggplot(data = estimates %>%
         filter(estimand %in% c("mu_1", "mu_0")) %>%
         filter(theta %in% thetas) %>%
         mutate(theta = paste0("alpha == ", theta)),
       aes(x = tau, y = est, group = estimand, color = estimand)) +
  geom_line(aes(linetype = estimand)) +
  facet_grid(. ~ theta, labeller = label_parsed) +
  labs(x = "Time (Days)", y = "Risk of Cholera") +
  scale_linetype_manual(values = c("mu_0" = "solid", "mu_1" = "dashed"),
                        labels = c("mu_0" = "Unvaccinated", "mu_1" = "Vaccinated"),
                        guide = guide_legend(title = NULL)) +
  scale_color_manual(values = c("mu_0" = "black", "mu_1" = "blue"),
                     labels = c("mu_0" = "Unvaccinated", "mu_1" = "Vaccinated"),
                     guide = guide_legend(title = NULL)) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(
  filename = "QuickStart_risk_over_time.png",
  width = 8,
  height = 4,
  dpi = 300
)

## Effects over theta plot
times = c(180,360)

ggplot(estimates %>%
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
    estimand = as_labeller(
      c("mu"   = "mu['B'](alpha)",
        "mu_1" = "mu['B,1'](alpha)",
        "mu_0" = "mu['B,0'](alpha)",
        "de"   = 'DE["B"](tau:alpha)',
        "se_1" = 'SE["B,1"](tau:alpha,0.45)',
        "se_0" = 'SE["B,0"](tau:alpha,0.45)',
        "oe"   = 'OE["B"](tau:alpha,0.45)'),
      label_parsed
    ),
    tau = label_parsed
  )) +
  
  scale_fill_manual(name = NULL,
                    values = c("95% Pointwise Confidence Interval" = "lightblue3",
                               "95% Uniform Confidence Band" = "lightblue1"),
                    guide = guide_legend(title = NULL)) +
  scale_color_manual(name = NULL,
                     values = c("Estimate" = "black"),
                     guide = guide_legend(title = NULL)) +
  labs(x = expression(alpha),
       y = "Effects") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank())

ggsave(
  filename = "QuickStart_effects_over_theta.png",
  width = 8,
  height = 6,
  dpi = 300
)
