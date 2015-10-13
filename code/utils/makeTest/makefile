SRCS_f90d1 = \
subroutines.f90 \
modules.f90 \
main.f90

OBJS = \
subroutines.f90 \
modules.f90 \
main.f90

PROG = main
# FC = /export/apps/pgi10.9/linux86/10.9/bin/pgf90
# FC = /share/apps/openmpi-1.8.5/pgi-2015/bin/mpif90
FC = /usr/bin/mpif90
IDIR =
LFLAGS =
LIBS =

#Set Optimization or Debug Flags
#     -O3   = Full Optimizations
#     -g    = Debug Mode
OPTMZ = -O3 -O2 

#Set Compiler Flags
#       -c       = Compile Only Don't Link (Required)
CFLAGS = -c

VPATH = $(SRC_DIR_f90d1):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90d1))

#Command Telling Make to Compile Program
$(PROG): $(MODULES) $(OBJS)
	$(FC) $(OBJS) $(MODULES) $(OPTMZ) -o $(PROG)

#Command Telling Make to Compile Modules
$(MODULES): %.o: %.f90
	$(FC) $(CFLAGS) $(OPTMZ) $< -o $@

#Command Telling Make to Compile Object Files
$(OBJS): %.o: %.f90 $(MODULES)
	$(FC) $(CFLAGS) $(OPTMZ) $< -o $@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)

# Dependencies of files
subroutines.o: \
    subroutines.f90 \
    modules.o
modules.o: \
    modules.f90
main.o: \
    main.f90 \
    modules.o \
    subroutines.o \
