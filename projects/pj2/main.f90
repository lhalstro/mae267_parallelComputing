! MAE 267
! PROJECT 2
! LOGAN HALSTROM
! 23 OCTOBER 2015


! DESCRIPTION:  Solve heat conduction equation for single block of steel.
! To compile: mpif90 -o main -O3 modules.f90 plot3D_module.f90 subroutines.f90 main.f90
    ! makes executable file 'main'
    ! 'rm *.mod' afterward to clean up unneeded compiled files
! To run: ./main or ./run.sh or sbatch run.sh on hpc1


PROGRAM heatTrans
!     USE CLOCK
    USE CONSTANTS
    USE subroutines
    USE plot3D_module

    IMPLICIT NONE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! INITIALIZE VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! BLOCKS
    TYPE(BLKTYPE), ALLOCATABLE :: blocks(:)
    ! GRID
    TYPE(MESHTYPE) :: mesh
    ! ITERATION PARAMETERS
    ! Minimum Residual
    REAL(KIND=8) :: min_res = 0.00001D0
    ! Maximum number of iterations
    INTEGER :: max_iter = 1000000, iter = 0, IBLK

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
    ! INIITIALIZE SOLUTION
    WRITE(*,*) 'Making mesh...'
    CALL init(blocks, mesh)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SOLVER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    WRITE(*,*) 'Solving heat conduction...'
!     CALL solve(mesh, min_res, max_iter, iter)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SAVE RESULTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    WRITE(*,*) 'Writing results...'
    ! SAVE SOLUTION AS PLOT3D FILES
    CALL plot3D(blocks)
    ! CALC TOTAL WALL TIME
    end_total = MPI_Wtime()
    wall_time_total = end_total - start_total
    ! SAVE SOLVER PERFORMANCE PARAMETERS
!     CALL output(mesh, iter)

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

    ! MOVE OUTPUT FILE TO OUTPUT DIRECTORY
    CALL EXECUTE_COMMAND_LINE ("mv a.out " // casedir // '.')

    WRITE(*,*) 'Done!'


END PROGRAM heatTrans
