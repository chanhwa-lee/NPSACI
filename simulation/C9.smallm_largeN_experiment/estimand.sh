#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8       # specify number of CPUs to use here
#SBATCH --array=1-100           # Number of Rscript submitting

module add r/4.5.0

# Parameters
M=100 # Number of clusters in one simulation
r=100  # Number of binary vector sampling
N=$1   # Cluster size (fixed)

# Define and create directories for output
mkdir -p "estimand/N${N}"

# Redirect stdout and stderr to custom log file
logfile="estimand/N${N}/estimand_rep-${SLURM_ARRAY_TASK_ID}.out"
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
  -N $N

# Add separator and job end message
echo "###################################################"
echo "### Job finished at $(date) ###"
echo "###################################################"