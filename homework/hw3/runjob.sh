#!/bin/bash -l
#SBATCH -J r_adia
#SBATCH -o r_adia-%J.out
#SBATCH -e r_adia-%J.err
#SBATCH -n number_processors

echo "starting at `date` on `hostname`"

#Run Command (assumes pgi and openmpi)

module load pgi openmpi hwloc
mpirun -n number_processors ./executable_name > "listing"
