! MAE 267
! PROJECT 3
! LOGAN HALSTROM
! 03 NOVEMBER 2015

! DESCRIPTION:  Subroutines used for solving heat conduction of steel plate.
! Subroutines utilizing linked lists are here so that linked lists do not need
! to be function inputs.
! Utilizes modules from 'modules.f90'
! CONTENTS:
! init --> Initialize the solution with dirichlet B.C.s
! solve --> Solve heat conduction equation with finite volume scheme
! output --> Save solution parameters to file

MODULE subroutines
    USE CONSTANTS
    USE BLOCKMOD
    USE IO

    IMPLICIT NONE

CONTAINS
    SUBROUTINE init_gridsystem(blocks)
        ! Initialize the solution with dirichlet B.C.s.  Save to restart files.
        TYPE(BLKTYPE)  :: blocks(:)

        ! INITIALIZE BLOCKS
        CALL init_blocks(blocks)
        ! WRITE BLOCK CONNECTIVITY FILE
        CALL write_blocks(blocks)
        ! INITIALIZE MESH
        CALL init_mesh(blocks)
        ! INITIALIZE TEMPERATURE WITH DIRICHLET B.C.
        CALL init_temp(blocks)

    END SUBROUTINE init_gridsystem

    SUBROUTINE init_solution(blocks, nbrlists)
        ! Read initial conditions from restart files.  Then calculate parameters
        ! used in solution
        TYPE(BLKTYPE)  :: blocks(:)
        ! LINKED LISTS STORING NEIGHBOR INFO
        TYPE(NBRLIST) :: nbrlists

        ! READ GRID


        ! READ INITIAL TEMPERATURE


        ! CALC LOCAL BOUNDARIES OF CELLS
        write(*,*) 'set local bounds'
        CALL set_block_bounds(blocks)
        ! INITIALIZE LINKED LISTS CONTAINING BOUNDARY INFORMATION
        write(*,*) 'make linked lists'
        CALL init_linklists(blocks, nbrlists)
        ! POPULATE BLOCK GHOST NODES
        write(*,*) 'update ghosts'
        CALL update_ghosts(blocks, nbrlists)
        ! CALC AREAS FOR SECONDARY FLUXES
        write(*,*) 'calc solution stuff'
        CALL calc_cell_params(blocks)
        ! CALC CONSTANTS OF INTEGRATION
        CALL calc_constants(blocks)

    END SUBROUTINE init_solution


!     SUBROUTINE solve(blocks, min_res, max_iter, iter)
!         ! Solve heat conduction equation with finite volume scheme
!         TYPE(BLKTYPE) :: blocks(:)
!         ! Minimum residual criteria for iteration, actual residual
!         REAL(KIND=8) :: min_res, res = 1000.D0
!         ! iteration number, maximum number of iterations
!         ! iter in function inputs so it can be returned to main
!         INTEGER :: iter, max_iter

!         INCLUDE "mpif.h"
!         REAL(KIND=8) :: start_solve, end_solve
!         WRITE(*,*) 'Starting clock for solver...'
!         start_solve = MPI_Wtime()

!         iter_loop: DO WHILE (res >= min_res .AND. iter <= max_iter)
!             ! Iterate FV solver until residual becomes less than cutoff or
!             ! iteration count reaches given maximum

!             ! INCREMENT ITERATION COUNT
!             iter = iter + 1
!             ! CALC NEW TEMPERATURE AT ALL POINTS
!             CALL calc_temp(blocks)

!             ! CALC RESIDUAL
!             res = MAXVAL( ABS( mesh%Ttmp(2:IMAX-1, 2:JMAX-1) ) )
!         END DO iter_loop

!         ! CACL SOLVER WALL CLOCK TIME
!         end_solve = MPI_Wtime()
!         wall_time_solve = end_solve - start_solve

!         ! SUMMARIZE OUTPUT
!         IF (iter > max_iter) THEN
!           WRITE(*,*) 'DID NOT CONVERGE (NUMBER OF ITERATIONS:', iter, ')'
!         ELSE
!           WRITE(*,*) 'CONVERGED (NUMBER OF ITERATIONS:', iter, ')'
!         END IF
!     END SUBROUTINE solve

    SUBROUTINE output(mesh, iter)
        ! Save solution parameters to file
        TYPE(MESHTYPE), TARGET :: mesh
        REAL(KIND=8), POINTER :: Temperature(:,:), tempTemperature(:,:)
        INTEGER :: iter, i, j

        Temperature => mesh%T(2:IMAX-1, 2:JMAX-1)
        tempTemperature => mesh%Ttmp(2:IMAX-1, 2:JMAX-1)
        ! Write final maximum residual and location of max residual
        OPEN(UNIT = 1, FILE = casedir // "SteadySoln.dat")
        DO i = 1, IMAX
            DO j = 1, JMAX
                WRITE(1,'(F10.7, 5X, F10.7, 5X, F10.7, I5, F10.7)'), mesh%x(i,j), mesh%y(i,j), mesh%T(i,j)
            END DO
        END DO
        CLOSE (1)

        ! Screen output
        WRITE (*,*), "IMAX/JMAX", IMAX, JMAX
        WRITE (*,*), "iters", iter
        WRITE (*,*), "residual", MAXVAL(tempTemperature)
        WRITE (*,*), "ij", MAXLOC(tempTemperature)

        ! Write to file
        OPEN (UNIT = 2, FILE = "SolnInfo.dat")
        WRITE (2,*), "Running a", IMAX, "by", JMAX, "grid took:"
        WRITE (2,*), iter, "iterations"
        WRITE (2,*), wall_time_total, "seconds (Total CPU walltime)"
        WRITE (2,*), wall_time_solve, "seconds (Solver CPU walltime)"
!         WRITE (2,*), wall_time_iter, "seconds (Iteration CPU walltime)"
        WRITE (2,*)
        WRITE (2,*), "Found max residual of ", MAXVAL(tempTemperature)
        WRITE (2,*), "At ij of ", MAXLOC(tempTemperature)
        CLOSE (2)
    END SUBROUTINE output



END MODULE subroutines


