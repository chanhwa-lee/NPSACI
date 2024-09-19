#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --output=log/estimator_rep-%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 # specify number of CPUs to use here

module load r/4.1.0

r=$1 # Number of sample splitting repetition
n=$2 # n_max: maximum cluster size
p=$3 # para or nonpara estimator

Rscript "../../estimator.R" \
  -r $r \
  -n $n \
  -p $p