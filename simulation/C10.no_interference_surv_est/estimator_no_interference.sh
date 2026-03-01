#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-1000
#SBATCH --export=NONE        # This prevents RStudio variables from leaking into the job

# 1. Start with a totally clean slate
module purge
unset R_HOME
unset R_LIBS

# 2. Load the RHEL9 version of R 4.4.0 (matching your RStudio)
module load r/4.5.0

# Parameters
m=$1 # Number of clusters per each simulation

# Define and create directories
d=estimate/No_interference_m${m}
mkdir -p $d/log $d/Rdata
cd $d

# Redirect logs
logfile="log/estimator_rep-${SLURM_ARRAY_TASK_ID}-${SLURM_JOB_ID}.out"
exec > >(tee -a "$logfile")
exec 2>&1


echo "##############################################################"
echo "### Job started at $(date) with m=${m} using no intereference methods ###"
echo "##############################################################"

# Run the R script with specified parameters
Rscript "../../estimator_no_interference.R" \
  -m $m 
  
# Add separator and job end message
echo "\n###################################################"
echo "### Job finished at $(date) ###"
echo "###################################################"