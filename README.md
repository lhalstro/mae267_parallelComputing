# Parallel Computations in Fluid/Thermal Sciences

##### MAE 267 - Fall Quarter 2015<br>UC Davis<br>Prof: Roger Davis

Parallel computing methods for solving 2D heat transfer with explicity, finite-volume schemes using Fortran and MPI.

### Projects
1. Serial, Single-Block Heat Conduction Solver
  * Solve heat conduction numerically on a single, non-uniform, 2D grid
  * Use one processor
2. Serial, Multi-Block Solver Intitialization
  * Initialize solver solution for single processor with grid decomposed into sub-domains (blocks)
  * Write block configuration files to restart solution
3. Serial, Multi-Block Solver
  * Modify solver to solve heat conduction on multiple blocks
  * Use ghost nodes to store neighbor information
4. Parallel, Multi-Block Solver Initialization
  * Initialize solver solution for multiple processors with grid decomposed into sub-domains
  * Distribute blocks to processors to achive equal load balancing
  * Save configuration files for each processor to restart solution
5. Parallel, Multi-Block Solver
  * Modify solver to solve heat conduction on multiple blocks with multiple processors
  * Include MPI code for multiple processors

### How To:

###### TO COMPILE A FORTRAN CODE (.f90 extension)

    GNU Fortran:
    'gfortran -o outname progname.f90'

    Open MPI (parallel, optimized):
    'mpif90 -o outname -O3 progname.f90'

    Using project specific makefile:
    './make.sh'

###### INSTALL OPEN-MPI WITH FORTRAN SUPPORT
    Homebrew tends to install open-mpi without fortran support.  To Fix:

    'brew reinstall openmpi --build-from-source'

###### RUNNING COMPILED CODES ON FRONT-END AND HPC1 BATCH SUBMISSION:
    
    MPI Parallel Run:
    'mpirun -n NCPUS main > "a.out"'
    
    Batch Job Submission on HPC1:
    'sbatch run.sh'

###### BASH ALIASES FOR SLURM QUEUING SYSTEM:
    alias fortmain='gfortran -o main main.f90'
    alias q="squeue"
    alias qq="squeue -u lhalstro"
    alias qdel="scancel"
    alias qsub="sbatch"

###### TO COPY FILES FROM REMOTE TO LOCAL (in local terminal):

    'scp username@remoteaddress:path/to/file path/to/copy/to'

    use -r to copy folders, i.e.:
    'scp -r lhalstro@hpc1.cse.ucdavis.edu:mae267/projects/pj5/Results projects/ucd/mae267/projects/pj5/.'

###### TO CHANGE FROM LOWER TO UPPER CASE IN SUBLIME TEXT:

    'ctrl+k,u' (l for lower)
