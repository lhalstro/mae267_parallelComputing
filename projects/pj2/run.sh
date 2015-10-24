#!/bin/bash -l
#SBATCH -J PJ1MAE267
#SBATCH -o slurm-%J.out
#SBATCH -e slurm-%J.err
NPROCS=1
#SBATCH -n $NPROCS

echo "starting at `date` on `hostname`"

#'&' will run process in background until it is complete
./main > "a.out"
