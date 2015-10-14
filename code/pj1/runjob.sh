#!/bin/bash -l
#SBATCH -J r_adia
#SBATCH -o r_adia-%J.out
#SBATCH -e r_adia-%J.err
#SBATCH -n 1

echo "starting at `date` on `hostname`"

#Run Command (assumes pgi and openmpi)

module load pgi openmpi hwloc
./main > "listing"
# mpirun -n number_processors ./main > "listing"
