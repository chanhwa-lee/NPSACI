# NPSACI

- The code in this repository implements the nonparametric cross-fitting estimators developed in the paper "**Nonparametric Causal Survival Analysis with Clustered Interference**"

## Quick start

Run `QuickStart.R` to quick start.

First, prepare dataset with column names 
- id: cluster id
- Y: observed time
- D: event indicator (D = I(T \le C))
- A: binary treatment status
- X: can be artibrary variable names indicating covariates

```r
# --- Generate toy example data ---
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
```

Next, specify 

```r
## Help functions for estimator main functions
source("/code/help_util.R")

## Help functions for policy specific functions
source("/code/help_TypeB.R")
source("/code/help_TPB.R")

## Help functions for nuisance functions estimation method
source("/code/help_nuis_est.R")

## Compute estimates
result = estimator(data = toy_data,
                   X.T.names = c("age", "dist.river"),
                   X.C.names = c("age"),
                   X.A.names = c("age", "dist.river"),
                   taus = 10*(1:50), thetas = seq(0.3, 0.6, length.out = 121), theta0 = 0.45)
```
 
## Summary

The GitHub repository comprises four folders: :file_folder:`code`, :file_folder:`simulation`, :file_folder:`application`, and :file_folder:`application_example`, which store reproducible codes for all computational work in the manuscript.

## :file_folder: code

Main R functions to implement proposed estimators.

- :page_facing_up: `help_util.R`: R functions for estimand and estimator computation
- :page_facing_up: `help_TypeB.R`: R functions specific to TypeB policy
- :page_facing_up: `help_TPB.R`: R functions specific to TPB policy
- :page_facing_up: `help_nuis_est.R`: R functions to fit and evaluate nuisance function estimators

## :file_folder: simulation

R files are the main scripts, while bash files (.sh) are for submitting parellel jobs to high-performance computing clusters (HPC) using SLURM.

- :page_facing_up: `help_simul.R`: Specifies the simulation setting (data-generating process)
- :page_facing_up: `estimand.R`: Computes target causal estimands
- :page_facing_up: `estimator.R`: Computes proposed estimators
- :page_facing_up: `readresult.R`: Read and summarize estimation result
- :file_folder: `estimand`: Target causal estimands
- :file_folder: `estimate`: Estimates of target causal estimands

## :file_folder: application

Analysis on the cholera vaccine effect on time to cholera incident accounting for the right censoring and clustered interference.
The raw data and computed estimates are not available to public.

- :page_facing_up: `preprocessing.R`: Preprocesses the raw data file and perform explarotry analysis.
- :page_facing_up: `estimator.R`: Computes proposed estimators
- :page_facing_up: `visualization.R`: Read and summarize estimation results to generate Figures


## :file_folder: application_example

Provides toy example to perform real data analysis instead of cholera vaccine data

- :page_facing_up: `generate_toy_data.R`: Generates toy data set mimicking the cholera vaccine data
- :file_folder: `estimate`: Estimates of target causal estimands using toy example dataset

***




# Detailed File Description

## :file_folder: simulation

- :file_folder: `M.main_simulation`: Perform simulations in main text and supplementary material C.1--C.4
- :file_folder: `C5.sigmab_experiment`: Perform simulations in supplementary material C.5
- :file_folder: `C6.Ndist_experiment`: Perform simulations in supplementary material C.6



R files are the main scripts, while bash files (.sh) are for submitting jobs to computing clusters.

### :file_folder: CIPS

CIPS policy simulation code and results included in the main text

- :page_facing_up: `Helpfunc.R`: R functions for estimand and estimator computation
- :page_facing_up: `estimand.R`: Causal estimands approximation
- :page_facing_up: `estimator.R`: Proposed estimators computation
- :page_facing_up: `readresult.R`: Read and summarize simulation result, reproduces Table 1
- :file_folder: `estimand`: Target causal estimands
- :file_folder: `data`: Estimates of target causal estimands

### :file_folder: TPB

TPB policy simulation code and results included in the main text

- :page_facing_up: `Helpfunc.R`: R functions for estimand and estimator computation
- :page_facing_up: `estimand.R`: Causal estimands approximation
- :page_facing_up: `estimator.R`: Proposed estimators computation
- :page_facing_up: `readresult.R`: Read and summarize simulation result, reproduces Table 2
- :page_facing_up: `estimands.RDS`: Target causal estimands
- :file_folder: `data`: Estimates of target causal estimands
  
### :file_folder: additional_simulation

Additional simulation code and results included in the supplementary material

> #### :file_folder: C.1.sizevarydelta
> CIPS with varying delta policy simulation code and results
> - :page_facing_up: `Helpfunc.R`: R functions for estimand and estimator computation
> - :page_facing_up: `estimand.R`: Causal estimands approximation
> - :page_facing_up: `estimator.R`: Proposed estimators computation
> - :page_facing_up: `readresult.R`: Read and summarize simulation result, reproduces Table S2
> - :file_folder: `estimand`: Target causal estimands
> - :file_folder: `data`: Estimates of target causal estimands

> #### :file_folder: C.2.Pro_IPW_comparison
> Proposed nonparametric estimator versus [Barkley et al. (2020)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-14/issue-3/Causal-inference-from-observational-studies-with-clustered-interference-with-application/10.1214/19-AOAS1314.full)'s IPW estimator comparison
> - :page_facing_up: `Helpfunc.R`: R functions for estimand and estimator computation
> - :page_facing_up: `compute_alpha.R`: Compute alpha values in Barkley et al. (2020) corresponding to delta values in CIPS
> - :page_facing_up: `alphas.rds`: Computed alpha values
> - :page_facing_up: `estimator.R`: Proposed estimators computation
> - :page_facing_up: `readresult.R`: Read and summarize simulation result, reproduces Figure S1
> - :file_folder: `data`: Estimates of target causal estimands 
> - NOTE: estimands were not computed here and instead loaded from `~/simulation/CIPS/estimand`

> #### :file_folder: C.3.r_comparison
> Proposed estimator performance over r (subsampling approximation degree) values
> - :page_facing_up: `Helpfunc.R`: R functions for estimand and estimator computation
> - :page_facing_up: `estimator.R`: Proposed estimators computation
> - :page_facing_up: `readresult.R`: Read and summarize simulation result, reproduces Figure S2
> - :file_folder: `data`: Estimates of target causal estimands
> - NOTE: estimands were not computed here and instead loaded from `~/simulation/CIPS/estimand`

> #### :file_folder: C.4.N_dist_comparison
> Proposed estimator performance over various cluster sizes distributions
> - :file_folder: `N3`, `N3_5`, `N5`, `N5_10`: Simulation results stored for various cluster sizes distributions. Structure is similar to `~/simulation/CIPS`
> - :page_facing_up: `readresult.R`: Read and summarize simulation result, reproduces Figure S3





## :file_folder: application

Analysis on the effect of water, sanitation, and hygiene (WASH) facilities on diarrhea among children in Senegal accounting for the clustered interference.

### :file_folder: Data

#### Senegal DHS data [(ANSD and ICF, 2020)](https://www.dhsprogram.com/pubs/pdf/FR368/FR368.pdf)
- Sociodemographic, enviromental, and health-related survey on household members 
- Used to assess the effect of WASH facilities on diarrhea incidence among children, allowing for interference within census blocks
- The data can be downloaded from the Demographic and Health Surveys Program website [https://dhsprogram.com](https://dhsprogram.com) after submitting a data request for research purposes. 
The process is as follows:

- **Register**
First, register for an account at the following link: [https://dhsprogram.com/data/new-user-registration.cfm](https://dhsprogram.com/data/new-user-registration.cfm). 
Fill in user information and select "Sub-Saharan Africa" from the drop-down menu under "Select Region". 
Then, click on the "Survey" checkbox for "Senegal" and submit the dataset request.

- **Download dataset**
Upon approval of your account registration, log in to the website at [https://dhsprogram.com/data/dataset_admin/login_main.cfm](https://dhsprogram.com/data/dataset_admin/login_main.cfm).
Select your project, then click on the "Download by Single Survey" link and select the country: "Senegal". 
Click on the "Download" link under the "Survey Datasets" column for "Country/Year: Senegal 2015" and "Type: Continuous DHS" row. 
Download "SNKR7HDT.ZIP" (Stata dataset) under the "Household Recode" tab. 
Uncompress the downloaded .ZIP folder. 
Rename "SNKR7HFL.DTA" to "senegal15.DTA". 
In the same folder, the .MAP file is the data dictionary.

- For other years (2016 - 2019), repeat the process according to the following steps:

  >- Senegal: Continuous DHS, 2015 -> (download) `SNKR7HDT.ZIP` -> (uncompress) `SNKR7HFL.DTA` -> (rename) `senegal15.DTA`
  >- Senegal: Continuous DHS, 2016 -> (download) `SNKR7IDT.ZIP` -> (uncompress) `SNKR7IFL.DTA` -> (rename) `senegal16.DTA`
  >- Senegal: Continuous DHS, 2017 -> (download) `SNKR7ZDT.ZIP` -> (uncompress) `SNKR7ZFL.DTA` -> (rename) `senegal17.DTA`
  >- Senegal: Continuous DHS, 2018 -> (download) `SNKR81DT.ZIP` -> (uncompress) `SNKR81FL.DTA` -> (rename) `senegal18.DTA`
  >- Senegal: Continuous DHS, 2019 -> (download) `SNKR8BDT.ZIP` -> (uncompress) `SNKR8BFL.DTA` -> (rename) `senegal19.DTA`

- Place the datasets in `~/application/Data/`

- :page_facing_up: `Preprocessing.R`: Preprocessing raw data files to generate `HHData.Rdata` and generate exploratory figures, reproduces Figures S4 and S5.

### :file_folder: CIPS

CIPS policy application code and results

- :page_facing_up: `estimator.R`: Proposed estimators computation
- :page_facing_up: `Visualization.R`: Read and summarize simulation results to generate Figures 1, S6, S7
- :file_folder: `Rdata`: Estimates of target causal estimands were stored.
- :file_folder: D.4. Comparison with [Park et al (2021)](https://arxiv.org/abs/2111.09932v1)
  >- :page_facing_up: `Preprocessing.R`: Preprocessing raw data files to generate HHData.Rdata from `~/application/Data/senegal18.DTA`
  >- :page_facing_up: `estimator.R`: Proposed estimators computation
  >- :page_facing_up: `Visualization.R`: Read and summarize simulation results to generate Figure S9
  >- :file_folder: `Rdata`: Estimates of target causal estimands

### :file_folder: TPB

TPB policy application code and results

- :page_facing_up: `Estimator.R`: Proposed estimators computation
- :page_facing_up: `Estimation.R`: Script for job submission to computing clusters
- :page_facing_up: `Visualization.R`: Read and summarize simulation results to generate Figures 2 and S8.
- :file_folder: `result`: Estimates of target causal estimands

***

## References
- [Barkley, B. G., Hudgens, M. G., Clemens, J. D., Ali, M. & Emch, M. E. (2020), ‘Causal
inference from observational studies with clustered interference, with application to a
cholera vaccine study’, The Annals of Applied Statistics 14(3), 1432–1448.](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-14/issue-3/Causal-inference-from-observational-studies-with-clustered-interference-with-application/10.1214/19-AOAS1314.full)

- [Park, C., Chen, G., Yu, M. & Kang, H. (2021), ‘Optimal allocation of water and sanitation
facilities to prevent communicable diarrheal diseases in Senegal under partial interference’,
arXiv preprint arXiv:2004.08950 .](https://arxiv.org/abs/2111.09932v1)https://arxiv.org/abs/2111.09932v1

