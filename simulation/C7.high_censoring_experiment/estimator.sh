#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-1000
#SBATCH --export=NONE        # This prevents RStudio variables from leaking into the job

# 1. Start with a totally clean slate
module purge
unset R_HOME
unset R_LIBS

# 2. Load the RHEL9 version of R 4.4.0 (matching your RStudio)
module load r/4.4.0

# Parameters
m=$1 # Number of clusters per each simulation
r=$2 # Number of subsampling approximation
c=$3 # Target censoring rate in simulated data ()

# Define and create directories
d=estimate/TypeB_m${m}_r${r}_CensoringRate${c}
mkdir -p $d/log $d/Rdata
cd $d

# Redirect logs
logfile="log/estimator_rep-${SLURM_ARRAY_TASK_ID}-${SLURM_JOB_ID}.out"
exec > >(tee -a "$logfile")
exec 2>&1


echo "##############################################################"
echo "### Job started at $(date) with m=${m}, r=${r}, CensoringRate=${c} ###"
echo "##############################################################"

# Run the R script with specified parameters
Rscript "../../estimator.R" \
  -m $m \
  -r $r \
  -c $c
  
# Add separator and job end message
echo "\n###################################################"
echo "### Job finished at $(date) ###"
echo "###################################################"