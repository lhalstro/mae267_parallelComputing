! MAE 267
! PROJECT 4
! LOGAN HALSTROM
! 14 NOVEMBER 2015


! DESCRIPTION:  Solve heat conduction equation for single block of steel.
! To compile:
    ! mpif90 -o main -O3 modules.f90 inout.f90 subroutines.f90 main.f90
        ! makes executable file 'main'
        ! 'rm *.mod' afterward to clean up unneeded compiled files
! To run on hpc1 nodes: sbatch run.sh
! To run on hpc1 front end: ./main or ./run.sh



PROGRAM heatTrans
!     USE CLOCK
    USE CONSTANTS
    USE subroutines
    USE IO

    IMPLICIT NONE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! INITIALIZE VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! BLOCKS
    TYPE(BLKTYPE), ALLOCATABLE :: blocks(:)
    ! LINKED LISTS STORING NEIGHBOR INFO
    TYPE(NBRLIST) :: nbrlists
    ! ITERATION PARAMETERS
    ! Residual history linked list
    TYPE(RESLIST), POINTER :: res_hist
    ! Maximum number of iterations
    INTEGER :: iter = 1, IBLK

    INCLUDE "mpif.h"
    REAL(KIND=8) :: start_total, end_total
    REAL(KIND=8) :: start_solve, end_solve
    ! CLOCK TOTAL TIME OF RUN
    start_total = MPI_Wtime()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! INITIALIZE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! READ INPUTS FROM FILE
    CALL read_input()
    ALLOCATE( blocks(NBLK) )
    ! INIITIALIZE GRID SYSTEM
    WRITE(*,*) 'Making mesh...'
    CALL init_gridsystem(blocks)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SOLVER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! INITIALIZE SOLUTION
    CALL init_solution(blocks, nbrlists)
    ! SOLVE
    WRITE(*,*) 'Solving heat conduction...'
    CALL solve(blocks, nbrlists, iter, res_hist)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SAVE RESULTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    WRITE(*,*) 'Writing results...'
    ! SAVE SOLUTION AS PLOT3D FILES
    CALL plot3D(blocks)
    ! CALC TOTAL WALL TIME
    end_total = MPI_Wtime()
    wall_time_total = end_total - start_total

    ! SAVE RESIDUAL HISTORY
    CALL write_res(res_hist)
    ! SAVE SOLVER PERFORMANCE PARAMETERS
    CALL output(blocks, iter)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! CLEAN UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO IBLK = 1, NBLK
        DEALLOCATE( blocks(IBLK)%mesh%xp   )
        DEALLOCATE( blocks(IBLK)%mesh%yp   )
        DEALLOCATE( blocks(IBLK)%mesh%x    )
        DEALLOCATE( blocks(IBLK)%mesh%y    )
        DEALLOCATE( blocks(IBLK)%mesh%T    )
        DEALLOCATE( blocks(IBLK)%mesh%Ttmp )
        DEALLOCATE( blocks(IBLK)%mesh%dt   )
        DEALLOCATE( blocks(IBLK)%mesh%V  )
        DEALLOCATE( blocks(IBLK)%mesh%V2nd )
        DEALLOCATE( blocks(IBLK)%mesh%term )
        DEALLOCATE( blocks(IBLK)%mesh%yPP)
        DEALLOCATE( blocks(IBLK)%mesh%yNP)
        DEALLOCATE( blocks(IBLK)%mesh%yNN)
        DEALLOCATE( blocks(IBLK)%mesh%yPN)
        DEALLOCATE( blocks(IBLK)%mesh%xNN)
        DEALLOCATE( blocks(IBLK)%mesh%xPN)
        DEALLOCATE( blocks(IBLK)%mesh%xPP)
        DEALLOCATE( blocks(IBLK)%mesh%xNP)
    END DO

    WRITE(*,*) 'Done!'

    ! MOVE OUTPUT FILE TO OUTPUT DIRECTORY
!     CALL EXECUTE_COMMAND_LINE ("mv a.out " // casedir // '.')


END PROGRAM heatTrans
