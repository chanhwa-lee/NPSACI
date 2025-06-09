# NPSACI

- The code in this repository implements the nonparametric cross-fitting estimators developed in the paper "**Nonparametric Causal Survival Analysis with Clustered Interference**"
 
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
- :file_folder: `Rdata`: Estimates of target causal estimands were stored.


## :file_folder: application_example

Provides toy example to perform real data analysis instead of cholera vaccine data






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



In the `~/application` folder, code for real data analysis is stored. 
Senegal Demographic Health Survey (DHS) provides sociodemographic, enviromental, and health-related information on household members. 
In the survey, households were randomly sampled from census blocks. 
The survey data was used to assess the effect of water, sanitation, and hygiene (WASH) facilities on diarrhea incidence among children in Senegal.
Senegal DHS WASH data is publicly available upon request at [https://dhsprogram.com/data/available-datasets.cfm](https://dhsprogram.com/data/available-datasets.cfm).
Place the Senegal DHS WASH data (described below) at `~/application/Data`, 
and run `~/application/Data/Preprocessing.R` to get cleaned data for the analysis. 
Run the `~/application/CIPS/estimator_main.sh` shell scripts to get estimates of Senegal DHS WASH causal estimands under the CIPS policy. The main code is in `estimator.R` in the same directory. The estimation results are summarized and visualized by `Visualization.R`. Estimation results for the TPB policy are done similarly, and the codes and results are in `~/application/TPB`.

***

## File Description

## :file_folder: code
- :page_facing_up: `help_util.R`: R functions for estimand and estimator computation
- :page_facing_up: `help_TypeB.R`: R functions specific to TypeB policy
- :page_facing_up: `help_TPB.R`: R functions specific to TPB policy
- :page_facing_up: `help_nuis_est.R`: R functions to fit and evaluate nuisance function estimators

## :file_folder: simulation

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

