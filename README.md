# NPSACI

This repository contains the code accompanying the paper:

* Lee, C., Zeng, D., Emch, M., Clemens, J. D., & Hudgens, M. G. (2024).
  *[Nonparametric Causal Survival Analysis with Clustered Interference](https://arxiv.org/abs/2409.13190)*

We develop nonparametric cross-fitting estimators for causal survival effects under clustered interference and right censoring.

---

## 📦 Quick Start

Run `QuickStart.R` to replicate a minimal example.

### Step 1: Prepare the dataset

Your dataset should contain the following variables:

* `id`: cluster ID
* `Y`: observed time
* `D`: event indicator (`D = I(T ≤ C)`)
* `A`: binary treatment
* `X`: covariates (any column names)

Here is a synthetic data generator:

```r
library(dplyr)
library(glue)
library(conflicted)
conflicts_prefer(dplyr::filter)

generate_cluster_data = function(N) {
  age <- round(runif(N, 15, 65))
  dist.river <- runif(N, 0, 5)
  b <- rnorm(1, 0, 0.5)
  pi <- plogis(0.2 + 0.2*(age/40 - 1)^2 + 0.2*pmax(dist.river/5, 0.3) + b)
  A <- rbinom(N, 1, pi)
  g.A <- (sum(A) - A) / (N - 1)
  T_ <- round(100 * rexp(N, 1 / exp(2 + A + g.A + A * g.A + 0.5 * dist.river - 0.5 * (age/40 - 1))))
  C <- round(runif(N, 0, 500 + 200 * A + 50 * dist.river + 100 * (age/40 - 1)))
  Y <- pmin(T_, C)
  D <- as.numeric(T_ <= C)
  data.frame(Y = Y, D = D, A = A, age = age, dist.river = dist.river)
}

m <- 100
toy_data <- lapply(sample(3:5, m, replace = TRUE), generate_cluster_data) %>%
  bind_rows(.id = "id") %>% mutate(id = as.numeric(id))
```

### Step 2: Estimate causal effects

Set parameters for nuisance models and target estimands:

* `X.T.names`: covariates used in the survival time model (`T`)
* `X.C.names`: covariates used in the censoring time model (`C`)
* `X.A.names`: covariates used in the treatment assignment model (`A`)
* `policy`: treatment allocation policy, either `"TypeB"` or `"TPB"`
* `taus`: time points of interest for survival analysis
* `thetas`: policy indices for which causal effects are estimated
* `theta0`: baseline policy index for causal effect comparison


```r
# --- Compute estimates ---

## Help functions for estimator main functions
source("code/help_util.R")

## Help functions for policy specific functions
source("code/help_TypeB.R")
source("code/help_TPB.R")

## Help functions for nuisance functions estimation method
source("code/help_nuis_est.R")

## Compute estimates
result <- estimator(
  data = toy_data,
  X.T.names = c("age", "dist.river"),
  X.C.names = c("age"),
  X.A.names = c("age", "dist.river"),
  policy = "TypeB",
  taus = 20 * (1:25),
  thetas = seq(0.3, 0.6, length.out = 31),
  theta0 = 0.45
)

## Estimation Result
result$result %>% filter(estimand == "mu", tau == 360) 
      estimand theta tau        est         se        PCL        PCU        UCL       UCU
#> 1        mu  0.30 360 0.17202130 0.04597763 0.08190515 0.26213745 0.07347707 0.2705655
#> 2        mu  0.31 360 0.16870090 0.04469703 0.08109472 0.25630708 0.07290139 0.2645004
#> 3        mu  0.32 360 0.16534002 0.04343936 0.08019888 0.25048117 0.07223608 0.2584440
#> 4        mu  0.33 360 0.16193982 0.04220139 0.07922510 0.24465455 0.07148923 0.2523904
#> 5        mu  0.34 360 0.15850180 0.04098049 0.07818003 0.23882356 0.07066796 0.2463356
#> 6        mu  0.35 360 0.15502781 0.03977458 0.07706964 0.23298598 0.06977862 0.2402770
#> 7        mu  0.36 360 0.15152010 0.03858205 0.07589929 0.22714091 0.06882688 0.2342133
#> 8        mu  0.37 360 0.14798130 0.03740175 0.07467387 0.22128873 0.06781781 0.2281448
#> 9        mu  0.38 360 0.14441442 0.03623297 0.07339781 0.21543104 0.06675600 0.2220728
#> 10       mu  0.39 360 0.14082292 0.03507536 0.07207521 0.20957062 0.06564560 0.2160002
```

**Result Columns**:

* `est`: estimated causal effect
* `se`: estimated standard error
* `PCL`/`PCU`: 95% **pointwise** confidence interval (L:lower, U:upper)
* `UCL`/`UCU`: 95% **uniform** confidence band (L:lower, U:upper)

---

## 📁 Repository Structure

### `/code`

Core functions for estimation:

* `help_util.R`: main estimation workflow
* `help_TypeB.R`: policy-specific functions (TypeB)
* `help_TPB.R`: policy-specific functions (TPB)
* `help_nuis_est.R`: nuisance parameter estimation

### `/simulation`

Scripts (.R) and HPC SLURM files (.sh) for simulation studies:

* `help_simul.R`: simulation settings
* `estimand.R`: computes target estimands
* `estimator.R`: computes estimators
* `readresult.R`: reads and summarizes outputs

**Subfolders:**

* `estimand/`: saved estimands
* `estimate/`: saved estimates

**Main Simulations in paper**

* **Performance under TypeB policy** (Table 1)
  
  :page_facing_up: `M.main_simulation/estimator.R`
  
  *Settings*: `policy = "TypeB"`, `m = 200`, `r = 100`

* **Proposed vs. IPCW estimators** (Figure 1)
  
  :page_facing_up: `M.main_simulation/estimator_chakladar.R`
  
  Compares the proposed estimator to the IPCW estimator by [Chakladar et al. (2021)](https://doi.org/10.1111/biom.13459)

**Supplementary Simulations**

* **Section C.1** – Performance under TPB policy (Table S1)
  
  :page_facing_up: `M.main_simulation/estimator.R`

  *Settings*: `policy = "TPB"`, `m = 200`, `r = 100`

* **Section C.2** – Performance over number of clusters *m* (Figure S1)

  :page_facing_up: `M.main_simulation/estimator.R`

  *Settings*: `policy = "TypeB"`, `m ∈ {25, 50, 100, 200, 400}`, `r = 100`

* **Section C.3** – Bounded vs. unbounded estimators (Figure S2)

  :page_facing_up: `M.main_simulation/estimator_unbounded.R`

  Evaluates estimators without the bounding modification

* **Section C.4** – Performance over subsampling degree *r* (Figure S3)

  :page_facing_up: `M.main_simulation/estimator.R`

  *Settings*: `policy = "TypeB"`, `m = 200`, `r ∈ {10, 20, 50, 100, 200, 500}`

* **Section C.5** – Performance under different correlation structures in treatment assignment (Figure S4)

  :file_folder: `C5.sigmab_experiment`

  Varies treatment correlation `σ_b = corr(A_ij, A_ik)`

* **Section C.6** – Performance under different cluster size distributions (Figure S5):

  :file_folder: `C6.Ndist_experiment`

  Varies distribution of cluster sizes `N_i`


### `/application`

Cholera vaccine effect analysis under clustered interference.
⚠️ **Raw data and outputs not public.**

* `preprocessing.R`: preprocessing + exploratory analysis (Figures S6–S8)
* `estimator.R`: causal estimation
* `visualization.R`: plotting (Figures 2–3, S9–S11)

### `/application_example`

Toy example replicating the application pipeline.

* `generate_toy_data.R`: create synthetic dataset
* `estimate/`: results using toy dataset

---

## 📚 References

* Lee, C., Zeng, D., Emch, M., Clemens, J. D., & Hudgens, M. G. (2024).
  *[Nonparametric Causal Survival Analysis with Clustered Interference](https://arxiv.org/abs/2409.13190)*. arXiv:2409.13190

* Chakladar, S. et al. (2022).
  *[Inverse probability weighted estimators of vaccine effects accommodating partial interference and censoring](https://doi.org/10.1111/biom.13459)*. Biometrics, 78(2), 777–788.
