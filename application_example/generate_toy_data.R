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

m = 100

toy_data = lapply(sample(x = 3:5, size = m, replace = T), generate_cluster_data) %>%
  dplyr::bind_rows(.id = "id") %>% mutate(id = as.numeric(id))

summary(toy_data)

write.csv(toy_data, "toy_data.csv", row.names = FALSE)

# --- Fit Cox proportional hazards model ---
model <- coxph(Surv(Y, D) ~ A + g.A + A*g.A + age + dist.river, 
               data = toy_data %>% group_by(id) %>% mutate(g.A = (sum(A)-A)/(n()-1)))

summary(model)