! MAE 267
! PROJECT 1
! LOGAN HALSTROM
! 12 OCTOBER 2015LS


! DESCRIPTION:  Solve heat conduction equation for single block of steel.
! To compile: mpif90 -o main modules.f90 plot3D_module.f90 subroutines.f90 main.f90

PROGRAM heatTrans
    USE CLOCK
    USE subroutines
    USE plot3D_module

    IMPLICIT NONE
    ! GRID
    TYPE(MESHTYPE), TARGET, ALLOCATABLE :: mesh(:,:)
    TYPE(CELLTYPE), TARGET, ALLOCATABLE :: cell(:,:)
    ! ITERATION PARAMETERS
    ! Minimum Residual
    REAL(KIND=8) :: min_res = 0.00001D0
    ! Maximum number of iterations
    INTEGER :: max_iter = 1000000, iter = 0

    ! MAKE GRID
    ! Set grid sizei
    CALL GRIDSIZE(101)
    ALLOCATE(mesh(1:IMAX, 1:JMAX))
    ALLOCATE(cell(1:IMAX-1, 1:JMAX-1))

    ! INIITIALIZE SOLUTION
    WRITE(*,*) 'Making mesh...'
    CALL init(mesh, cell)
    ! MEASURE WALL TIME FOR OVERALL SOLUTION
    WRITE(*,*) 'Starting clock for solver...'
    CALL start_clock()
    ! SOLVE
    WRITE(*,*) 'Solving heat conduction...'
    CALL solve(mesh, cell, min_res, max_iter, iter)
    CALL end_clock()
    ! SAVE SOLUTION PARAMETERS
    WRITE(*,*) 'Writing results...'
    CALL output(mesh, iter)
    ! SAVE SOLUTION AS PLOT3D FILES
    CALL plot3D(mesh)
    ! CLEAN UP
    DEALLOCATE(mesh)
    DEALLOCATE(cell)
    WRITE(*,*) 'Done!'

END PROGRAM heatTrans
