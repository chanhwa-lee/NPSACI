#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --output=estimand/estimand_rep-%a.out

M=1000 # Number of clusters in one simulation
r=100 # Number of binary vector sampling

Rscript "estimand.R" \
  -M $M \
  -r $r