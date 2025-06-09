###---- Generate synthetic data ----###

m = 100

generate_cluster_data = function(N){
  
  ## Step1. Covariate generation
  age <- round(runif(N, 15, 65))
  dist.river <- runif(N, 0, 5)
  
  ## Step2. Treatment model
  sigma.b = 0.5
  b = rnorm(1,0,sigma.b)
  pi = plogis(0.2 + 0.2*(age/40-1)^2 + 0.2*pmax(dist.river/5,0.3) + b)
  A <- rbinom(N, 1, pi)
  g.A <- (sum(A)-A)/(N-1)
  
  ## Step3. Survival time model
  
  ### Generate T from exponential distribution
  T_ <- round(100*rexp(N, 1/exp(2 + 1*A + 1*g.A + 1*A*g.A + 0.5*dist.river - 0.5*(age/40-1))))
  # T_ <- round(100*rexp(N, 1/exp(2 + 0.8*A + 1*g.A + 0.8*A*g.A + 0.5*dist.river - 0.5*(age/40-1))))
  
  ### Generate C from uniform distribution
  C <- round(runif(N, 0, 500 + 200*A + 50*dist.river + 100*(age/40-1)))
  
  ## Step 4. Combine all data
  Y <- pmin(T_,C)
  D <- as.numeric(T_ <= C)
  
  return(data.frame(Y = Y, T = T_, C = C, D = D, A = A, g.A = g.A, age = age, dist.river = dist.river))
  
}

toy_data = lapply(sample(x = 3:5, size = m, replace = T), generate_cluster_data)

toy_data = dplyr::bind_rows(toy_data,.id = "id") %>% mutate(id = as.numeric(id))

summary(toy_data)
summary(data)

# Fit Cox proportional hazards model
model <- coxph(Surv(Y, D) ~ A + g.A + age + dist.river, data = toy_data)

# View model summary
summary(model)

write.csv(toy_data, "toy_data.csv", row.names = FALSE)
