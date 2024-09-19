#!/bin/bash

module add r/4.1.0

mkdir -p estimand

### Estimand Simulation ###
jid_estimand=$( sbatch --array=1-50 estimand.sh | awk '{print $4}' )
