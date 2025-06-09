#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/dev/null
#SBATCH --array=1-100  # Number of Rscript submitting

module add r/4.1.0

# Parameters
M=1000 # Number of clusters in one simulation
r=100  # Number of binary vector sampling
p=$1   # Policy

# Define and create directories for output
mkdir -p "estimand/$p"

# Redirect stdout and stderr to custom log file
logfile="estimand/$p/estimand_rep-${SLURM_ARRAY_TASK_ID}.out"
exec > >(tee -a "$logfile")
exec 2>&1

# Add separator and job start message
echo "##############################################################"
echo "### Job started at $(date) ###"
echo "##############################################################"

# Run the R script with specified parameters
Rscript "estimand.R" \
  -M $M \
  -r $r \
  -p $p

# Add separator and job end message
echo "###################################################"
echo "### Job finished at $(date) ###"
echo "###################################################"