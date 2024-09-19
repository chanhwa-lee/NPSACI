#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=log/estimator_rep-%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 # specify number of CPUs to use here

module load r/4.1.0

m=$1 # Number of clusters per simulation
r=$2 # Number of sample splitting repetition
p=$3 # para or nonpara estimator
b=$4 # sigma.b

Rscript "../../estimator.R" \
  -m $m \
  -r $r \
  -p $p \
  -b $b