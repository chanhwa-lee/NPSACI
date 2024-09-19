#!/bin/bash

module add r/4.1.0

### Estimator Simulation ###
S=$1 # Number of sample splitting repetition
r=$2 # Number of subsampling approximation
n=$3 # n_max: maximum cluster size
p=$4 # para or nonpara estimator

d=estimate/S${S}_r${r}_nmax${n}_${p}
echo $d
mkdir -p $d
mkdir -p $d/log
mkdir -p $d/Rdata
cd $d

jid_estimator=$( sbatch --array=1-$S ../../estimator.sh $r $n $p | awk '{print $4}' )
