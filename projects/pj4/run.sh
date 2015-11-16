#!/bin/bash -l
#SBATCH -J PJ4MAE267
#SBATCH -o slurm-%J.out
#SBATCH -e slurm-%J.err
NPROCS=6
#SBATCH -n $NPROCS

#Text with start time and location
START="starting at `date` on `hostname`"
#Print to screen
echo $START
#Print to output file
echo $START > "a.out"

# '&' will run process in background until it is complete
mpirun -n $NPROCS main >> "a.out"
# ./main >> "a.out"

# Move files to case directory
python3 move.py
