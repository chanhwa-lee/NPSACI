#!/bin/bash
#SBATCH --time=05:30:00
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=8       # specify number of CPUs to use here
#SBATCH --array=1-100  # Number of Rscript submitting
#SBATCH --export=NONE        # This prevents RStudio variables from leaking into the job

# 1. Start with a totally clean slate
module purge
unset R_HOME
unset R_LIBS

module load r/4.5.0

# Parameters
m=$1 # Number of clusters per each simulation

# Define and create directories for output
d=nuis_estimates/m${m}
mkdir -p $d/log $d/Rdata
cd $d

# Redirect stdout and stderr to custom log file
logfile="log/estimator_rep-${SLURM_ARRAY_TASK_ID}.out"
exec > >(tee -a "$logfile")
exec 2>&1

# Add separator and job start message
echo "##############################################################"
echo "### Job started at $(date) with m=${m}  ###"
echo "##############################################################"

# Run the R script with specified parameters
Rscript "../../check_nuis_converge.R" \
  -m $m 

# Add separator and job end message
echo "###################################################"
echo "### Job finished at $(date) ###"
echo "###################################################"
