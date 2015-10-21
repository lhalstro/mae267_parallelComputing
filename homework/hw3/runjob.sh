#!/bin/bash -l
#SBATCH -J r_adia
#SBATCH -o r_adia-%J.out
#SBATCH -e r_adia-%J.err
#SBATCH -n 2s

echo "starting at `date` on `hostname`"

#Run Command (assumes pgi and openmpi)

module load pgi openmpi hwloc
mpirun -n 2 calcpip > "a.out"
