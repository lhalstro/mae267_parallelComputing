# MAE 267 - Parallel Computations in Fluid/Thermal Sciences

##### Fall Quarter 2015<br>UC Davis<br>Prof: Roger Davis

Paralell computing methods for solving 2D heat transfer with explicity, finite-volume schemes using Fortran and MPI.

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

    'gfortran -o outname progname.f90'

###### INSTALL OPEN-MPI WITH FORTRAN SUPPORT
    Homebrew tends to install open-mpi without fortran support.  To Fix:

    'brew reinstall openmpi --build-from-source'

###### TO COPY FILES FROM REMOTE TO LOCAL (in local terminal):

    'scp username@remoteaddress:path/to/file path/to/copy/to'

    use -r to copy folders, i.e.:
    'scp -r lhalstro@hpc1.cse.ucdavis.edu:mae267/code/pj1/results projects/ucd/mae267/code/pj1/.'

###### TO CHANGE FROM LOWER TO UPPER CASE IN SUBLIME TEXT:

    'ctrl+k,u' (l for lower)
