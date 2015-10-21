#!/bin/bash -l
#SBATCH -J calcpip
#SBATCH -o slurm-%J.out
#SBATCH -e slurm-%J.err
NPROCS=1
#SBATCH -n $NPROCS

echo "starting at `date` on `hostname`"

#Run Command (assumes pgi and openmpi)

#module load pgi openmpi hwloc
mpirun -n $NPROCS calcpip > "a.out"
