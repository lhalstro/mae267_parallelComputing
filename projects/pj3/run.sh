#!/bin/bash -l
#SBATCH -J PJ1MAE267
#SBATCH -o slurm-%J.out
#SBATCH -e slurm-%J.err
NPROCS=1
#SBATCH -n $NPROCS

#Text with start time and location
START="starting at `date` on `hostname`"
#Print to screen
echo $START
#Print to output file
echo $START >> "a.out"

# '&' will run process in background until it is complete
#./main >> "a.out"
./main
