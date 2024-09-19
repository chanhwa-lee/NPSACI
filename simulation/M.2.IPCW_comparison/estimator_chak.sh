#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=log/estimator_rep-%a.out
#SBATCH --ntasks=1

m=$1 # Number of clusters per simulation

module load r/4.1.0
Rscript "../../estimator_chak.R" -m $m