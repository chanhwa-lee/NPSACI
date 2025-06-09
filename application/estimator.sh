#!/bin/bash
#SBATCH --time=7-01:00:00
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 # specify number of CPUs to use here
#SBATCH --array=1-15  # Number of sample splitting repetition

module add r/4.1.0

# Parameters
n=$1 # n_max: maximum cluster size
r=$2 # Number of subsampling approximation
s=$3 # Scenario indicating which nuisance estimators are used
p=$4 # Policy

# Define and create directories for output
d=estimate/${p}_nmax${n}_r${r}_Scenario${s}
mkdir -p $d/log $d/Rdata
cd $d

# Redirect stdout and stderr to custom log file
logfile="log/estimator_rep-${SLURM_ARRAY_TASK_ID}-${SLURM_JOB_ID}.out"
exec > >(tee -a "$logfile")
exec 2>&1

# Add separator and job start message
echo "##############################################################"
echo "### Job started at $(date) with nmax=${n}, r=${r}, s=${s} ###"
echo "##############################################################"

# Run the R script with specified parameters
Rscript "../../estimator.R" \
  -n $n \
  -r $r \
  -s $s \
  -p $p

# Add separator and job end message
echo "###################################################"
echo "### Job finished at $(date) ###"
echo "###################################################"
