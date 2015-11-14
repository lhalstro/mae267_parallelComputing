! MAE 267
! PROJECT 3
! LOGAN HALSTROM
! 03 NOVEMBER 2015

! DESCRIPTION:  Subroutines used for solving heat conduction of steel plate.
! Subroutines utilizing linked lists are here so that linked lists do not need
! to be function inputs.
! Utilizes modules from 'modules.f90'

! CONTENTS:
    ! init_gridsystem
        ! Initialize the solution with dirichlet B.C.s.  Save to restart files.

    ! init_solution
        ! Read initial conditions from restart files.  Then calculate parameters
        ! used in solution

    ! solve
        ! Solve heat conduction equation with finite volume scheme
        ! (within iteration loop)

    ! output
        ! Save solution performance parameters to file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
        ! WRITE GRID AND INITIAL TEMPERATURE TO PLOT3D RESTART FILES
        CALL plot3D(blocks)

    END SUBROUTINE init_gridsystem

    SUBROUTINE init_solution(blocks, nbrlists)
        ! Read initial conditions from restart files.  Then calculate parameters
        ! used in solution

        TYPE(BLKTYPE)  :: blocks(:)
        ! LINKED LISTS STORING NEIGHBOR INFO
        TYPE(NBRLIST) :: nbrlists

        ! READ BLOCK CONFIGURATION INFORMATION FROM CONFIG FILE
        CALL read_blocks(blocks)

        ! READ GRID AND INITIAL TEMPERATURE FROM PLOT3D RESTART FILE
        CALL readPlot3D(blocks)


        ! CALC LOCAL BOUNDARIES OF CELLS
        write(*,*) 'set local bounds'
        CALL set_block_bounds(blocks)



        ! INITIALIZE LINKED LISTS CONTAINING BOUNDARY INFORMATION
        write(*,*) 'make linked lists'
        CALL init_linklists(blocks, nbrlists)
        ! POPULATE BLOCK GHOST NODES
        write(*,*) 'update ghosts'
        CALL update_ghosts(blocks, nbrlists)

!         CALL update_ghosts_debug(blocks)

        ! CALC AREAS FOR SECONDARY FLUXES
        write(*,*) 'calc solution stuff'
        CALL calc_cell_params(blocks)
        ! CALC CONSTANTS OF INTEGRATION
        CALL calc_constants(blocks)

    END SUBROUTINE init_solution


    SUBROUTINE solve(blocks, nbrlists, iter, res_hist)
        ! Solve heat conduction equation with finite volume scheme
        ! (within iteration loop)

        TYPE(BLKTYPE) :: blocks(:)
        ! LINKED LISTS STORING NEIGHBOR INFO
        TYPE(NBRLIST) :: nbrlists
        ! Residual history linked list
        TYPE(RESLIST), POINTER :: res_hist
        ! pointer to iterate linked list
        TYPE(RESLIST), POINTER :: hist
        ! Minimum residual criteria for iteration, actual residual
        REAL(KIND=8) :: res = 1000.D0, resloc, resmax
        ! iter in function inputs so it can be returned to main
        INTEGER :: iter, IBLK, IBLKRES

        INCLUDE "mpif.h"
        REAL(KIND=8) :: start_solve, end_solve
        WRITE(*,*) 'Starting clock for solver...'
        start_solve = MPI_Wtime()

        ! residual history
        ALLOCATE(res_hist)
        hist => res_hist

        iter_loop: DO WHILE (res >= min_res .AND. iter <= max_iter)
            ! Iterate FV solver until residual becomes less than cutoff or
            ! iteration count reaches given maximum

            ! CALC NEW TEMPERATURE AT ALL POINTS
            CALL calc_temp(blocks)

            ! UPDATE GHOST NODES WITH NEW TEMPERATURE SOLUTION
            CALL update_ghosts(blocks, nbrlists)
!             CALL update_ghosts_debug(blocks)

            ! CALC RESIDUAL
            resmax = 0.D0
            DO IBLK = 1, NBLK
                ! Find max of each block
                resloc = MAXVAL( ABS( blocks(IBLK)%mesh%Ttmp(2:IMAXBLK-1, 2:JMAXBLK-1) ) )
                ! keep biggest residual
                IF (resmax < resloc) THEN
                    resmax = resloc
                END IF
            END DO
            ! FINAL RESIDUAL
            res = resmax

            ! SWITCH TO NEXT LINK
                ! (skip first entry)
            ALLOCATE(hist%next)
            hist => hist%next
            NULLIFY(hist%next)
            ! STORE RESIDUAL HISTORY
            hist%iter = iter
            hist%res = res


            ! INCREMENT ITERATION COUNT
            iter = iter + 1

        END DO iter_loop

        ! there was an extra increment after final iteration we need to subtract
        iter = iter - 1

        ! CACL SOLVER WALL CLOCK TIME
        end_solve = MPI_Wtime()
        wall_time_solve = end_solve - start_solve

        IF (iter > max_iter) THEN
          WRITE(*,*) 'DID NOT CONVERGE (NUMBER OF ITERATIONS:', iter, ')'
        ELSE
          WRITE(*,*) 'CONVERGED (NUMBER OF ITERATIONS:', iter, ')'
          WRITE(*,*) '          (MAXIMUM RESIDUAL    :', res,  ')'
        END IF
    END SUBROUTINE solve

    SUBROUTINE output(blocks, iter)
        ! Save solution performance parameters to file

        TYPE(BLKTYPE), TARGET :: blocks(:)
        REAL(KIND=8), POINTER :: tmpT(:,:), tempTemperature(:,:)
        REAL(KIND=8) :: resloc, resmax
        INTEGER :: iter, I, J, IBLK, IRES

!         Temperature => mesh%T(2:IMAX-1, 2:JMAX-1)
!         tempTemperature => mesh%Ttmp(2:IMAX-1, 2:JMAX-1)

        ! CALC RESIDUAL
        resmax = 0.D0
        DO IBLK = 1, NBLK
            ! Find max of each block
            resloc = MAXVAL( ABS( blocks(IBLK)%mesh%Ttmp(2:IMAXBLK-1, 2:JMAXBLK-1) ) )
            ! keep biggest residual
            IF (resmax < resloc) THEN
                resmax = resloc
                IRES = IBLK
            END IF
        END DO


        ! Write final maximum residual and location of max residual
!         OPEN(UNIT = 1, FILE = casedir // "SteadySoln.dat")
!         DO i = 1, IMAX
!             DO j = 1, JMAX
!                 WRITE(1,'(F10.7, 5X, F10.7, 5X, F10.7, I5, F10.7)'), mesh%x(i,j), mesh%y(i,j), mesh%T(i,j)
!             END DO
!         END DO
!         CLOSE (1)

        ! Screen output
        tmpT => blocks(IRES)%mesh%Ttmp
        WRITE (*,*), "IMAX/JMAX", IMAX, JMAX
        WRITE (*,*), "N/M", N, M
        WRITE (*,*), "iters", iter
        WRITE (*,*), "max residual", MAXVAL(tmpT(2:IMAXBLK-1, 2:JMAXBLK-1))
        WRITE (*,*), "on block id", IRES
        WRITE (*,*), "residual ij", MAXLOC(tmpT(2:IMAXBLK-1, 2:JMAXBLK-1))

        ! Write to file
!         OPEN (UNIT = 2, FILE = TRIM(casedir) // "SolnInfo.dat")
        OPEN (UNIT = 2, FILE = "SolnInfo.dat")
        WRITE (2,*), "Running a", IMAX, "by", JMAX, "grid,"
        WRITE (2,*), "With NxM:", N, "x", M, "blocks took:"
        WRITE (2,*), iter, "iterations"
        WRITE (2,*), wall_time_total, "seconds (Total CPU walltime)"
        WRITE (2,*), wall_time_solve, "seconds (Solver CPU walltime)"
!         WRITE (2,*), wall_time_iter, "seconds (Iteration CPU walltime)"
        WRITE (2,*)
        WRITE (2,*), "Found max residual of ", MAXVAL(tmpT(2:IMAXBLK-1, 2:JMAXBLK-1))
        WRITE (2,*), "on block id", IRES
        WRITE (2,*), "At ij of ", MAXLOC(tmpT(2:IMAXBLK-1, 2:JMAXBLK-1))
        CLOSE (2)
    END SUBROUTINE output



END MODULE subroutines


