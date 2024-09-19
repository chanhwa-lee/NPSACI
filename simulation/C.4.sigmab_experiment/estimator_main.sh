#!/bin/bash

module add r/4.1.0

### Estimator Simulation ###
m=$1 # Number of clusters per each simulation
r=$2 # Number of subsampling approximation
p=$3 # para or nonpara estimator
b=$4 # sigma.b

d=estimate/m${m}_r${r}_${p}_sigmab${b}
echo $d
mkdir -p $d
mkdir -p $d/log
mkdir -p $d/Rdata
cd $d

jid_estimator=$( sbatch --array=1-1000 ../../estimator.sh $m $r $p $b | awk '{print $4}' )
