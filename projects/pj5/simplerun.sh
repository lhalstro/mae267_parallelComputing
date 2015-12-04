#!/bin/bash -l

#SBATCH -J lhalstro_PJ5_NP4
#SBATCH -o lhalstro_PJ5_NP4-%J.out
#SBATCH -e lhalstro_PJ5_NP4-%J.err
#SBATCH -n 4

echo "starting at `date` on `hostname`"

#Run Command (assumes pgi and openmpi)

module load openmpi hwloc
mpirun -n 4 ./main > "a.out"
