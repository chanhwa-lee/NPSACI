#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8       # specify number of CPUs to use here
#SBATCH --array=1-1000           # Number of Rscript submitting

module add r/4.1.0

# Parameters
m=$1 # Number of clusters per each simulation
r=$2 # Number of subsampling approximation
s=$3   # Correlation between Aij and Aik, sigma.b

# Define and create directories for output
d=estimate/TypeB_m${m}_r${r}_sigma.b${s}
mkdir -p $d/log $d/Rdata
cd $d

# Redirect stdout and stderr to custom log file
logfile="log/estimator_rep-${SLURM_ARRAY_TASK_ID}-${SLURM_JOB_ID}.out"
exec > >(tee -a "$logfile")
exec 2>&1

# Add separator and job start message
echo "##############################################################"
echo "### Job started at $(date) with m=${m}, r=${r}, sigma.b=${s} ###"
echo "##############################################################"

# Run the R script with specified parameters
Rscript "../../estimator.R" \
  -m $m \
  -r $r \
  -s $s

# Add separator and job end message
echo "###################################################"
echo "### Job finished at $(date) ###"
echo "###################################################"
