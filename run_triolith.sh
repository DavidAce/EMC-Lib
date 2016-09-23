#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-08:00:00
#SBATCH -N 8
#SBATCH --exclusive
export OMP_NUM_THREADS=1
mpprun ./build/Release/EMC
