MAE 267 - Parallel Computations in Fluid/Thermal Sciences
Prof. Roger Davis
Fall 2015

TO COMPILE A FORTRAN CODE (.f90 extension)
    'gfortran -o outname progname.f90'

INSTALL OPEN-MPI WITH FORTRAN SUPPORT
    Homebrew tends to install open-mpi without fortran support.  To Fix:
    'brew reinstall openmpi --build-from-source'

TO COPY FILES FROM REMOTE TO LOCAL (in local terminal):
    'scp username@remoteaddress:path/to/file path/to/copy/to'

    use -r to copy folders, i.e.:
    'scp -r lhalstro@hpc1.cse.ucdavis.edu:mae267/code/pj1/results projects/ucd/mae267/code/pj1/.'

TO CHANGE FROM LOWER TO UPPER CASE IN SUBLIME TEXT:
    'ctrl+k,u' (l for lower)
