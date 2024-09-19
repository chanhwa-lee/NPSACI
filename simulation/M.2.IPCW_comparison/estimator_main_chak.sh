#!/bin/bash

module add r/4.1.0

### Estimator Simulation ###
m=$1 # Number of clusters per each simulation
a=1000 # Simulation replicates

### Chakladar IPCW estimator ###
d=estimate_chak/m${m}
echo $d
mkdir -p $d
mkdir -p $d/log
mkdir -p $d/Rdata
cd $d

sbatch --array=1-$a ../../estimator_chak.sh $m
