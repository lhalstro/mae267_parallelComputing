# makefile for codename on RedHat Linux 6.2 workstations

# some utilities
ECHO    = echo
RM      = /bin/rm -f
CHMOD   = chmod
MPI_HOME= /opt/openmpi

# set source directory and MPI include directory
# souce directory currently set as current directory where
# makefile and code is located.  this could be changed to 
# allow source and objects (along with make) to be located in
# separate directories
codename_SRC = ./

# fortran-90 compiler and linker
# GNU compiler:
#FC    = mpif90
#FL    = mpif90
# Portland-Group compiler:
FC    = /share/apps/openmpi-pgi/bin/mpif90
FL    = /share/apps/openmpi-pgi/bin/mpif90

# C compiler
# GNU compiler:
#CC    = mpicc
# Portland-Group compiler:
CC    = /share/apps/openmip-pgi/bin/mpicc

# compiler flag for compiling only
PROC = -c 

# compiler flags for optimization/debug
FOPT =  
COPT =

# architecture on which fortran90 compiler are running
ARCH = 

# default fortran-90 compiler flags
FFLAGS = 

# default c compiler flags
CFLAGS = 

# form one flag variable to feed to compiler 
FFLAGS_codename = $(ARCH) $(FOPT) $(FFLAGS)
CFLAGS_codename = $(ARCH) $(COPT) $(CFLAGS)

# linker flag for executables 
PROL = -o

# fortran90 linker flags
LDFLAGS = $(ARCH)

# extra lib path and lib name
LIBS = 
LIBS_pV3 = 

# suffixes for this environment; object,executable,library,archive
O=.o
X=
M=.mod
L=.a
A=
AWKS=.awk
SEDS=.sed
LEXS=.l
YACS=.y

#======================================================================

# dependency file
MAKEFILE_PID = makefile.pid

# default targets used by all makefiles
default: build

# build codename 
build: codename$(X)

# Initialize these macros so they are defined.
ALL_OBJS_codename=
OBJ_PI_codename=
OBJ_PD_codename=
include $(MAKEFILE_PID)

ALL_OBJS_codename=$(OBJ_PI_codename) $(OBJ_PD_codename)

codename$(X): $(ALL_OBJS_codename)
	$(FL) $(PROL) codename$(X) $(LDFLAGS) $(ALL_OBJS_codename) \
 	$(CHMOD) 775 codename$(X)

# clean all modules and objects and binaries
clean:
	$(RM) codename$(X) *$(O) *$(M) core
