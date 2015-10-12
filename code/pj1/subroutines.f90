! MAE 267
! PROJECT 1
! LOGAN HALSTROM
! 12 OCTOBER 2015

! DESCRIPTION:  Subroutines used for solving heat conduction of steel plate.
! Utilizes modules from 'modules.f90'
! CONTENTS:
! init --> Initialize the solution with dirichlet B.C.s
! solve --> Solve heat conduction equation with finite volume scheme

MODULE subroutines
    USE CONSTANTS
    USE MAKEGRID
    USE CELLS
    USE TEMPERATURE

    IMPLICIT NONE

CONTAINS
    SUBROUTINE init(mesh, cells)
        ! Initialize the solution with dirichlet B.C.s
        TYPE(GRID), TARGET :: mesh(1:IMAX, 1:JMAX)
        TYPE(CELL), TARGET :: cells(1:IMAX-1, 1:JMAX-1)
        INTEGER :: i, j

        ! INITIALIZE MESH
        CALL init_mesh(mesh)
        ! INITIALIZE CELLS
        CALL init_cells(mesh, cells)
        ! CALC SECONDARY AREAS OF INTEGRATION
        CALL calc_2nd_areas(mesh, cells)
        ! CALC CONSTANTS OF INTEGRATION
        CALL calc_constants(mesh, cells)

        ! INITIALIZE TEMPERATURE WITH DIRICHLET B.C.
        !PUT DEBUG BC HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO j = 1, JMAX
          CALL init_temp(mesh(1,j), 3.D0 * mesh(1,j)%yp + 2.D0)
          CALL init_temp(mesh(IMAX,j), 3.D0 * mesh(IMAX,j)%yp + 2.D0)
        END DO

        DO i = 1, IMAX
          CALL init_temp(mesh(i,1), ABS(COS(pi * mesh(i,1)%xp)) + 1.D0)
          CALL init_temp(mesh(i,JMAX), 5.D0 * (SIN(pi * mesh(i,JMAX)%xp) + 1.D0))
        END DO
    END SUBROUTINE init

    SUBROUTINE solve(mesh, cells, min_res, max_iter)
        ! Solve heat conduction equation with finite volume scheme
        TYPE(GRID) :: mesh(1:IMAX, 1:JMAX)
        TYPE(CELL) :: cells(1:IMAX-1, 1:JMAX-1)
        ! Minimum residual criteria for iteration, actual residual
        REAL(KIND=8) :: min_res, res = 1000.D0
        ! iteration number, maximum number of iterations
        INTEGER :: iter = 0, max_iter
        INTEGER :: i, j

        iter_loop: DO WHILE (res >= min_res .AND. iter <= max_iter)
            ! Iterate FV solver until residual becomes less than cutoff or
            ! iteration count reaches given maximum

            ! INCREMENT ITERATION COUNT
            iter = iter + 1
            ! CALC NEW TEMPERATURE AT ALL POINTS
            CALL derivatives(mesh, cells)
            ! SAVE NEW TEMPERATURE DISTRIBUTION
            DO j = 2, JMAX - 1
                DO i = 2, IMAX - 1
                    mesh(i,j)%T = mesh(i,j)%T + mesh(i,j)%Ttmp
                END DO
            END DO
            ! CALC RESIDUAL
            res = MAXVAL(ABS(mesh(2:IMAX-1, 2:JMAX-1)%Ttmp))
        END DO iter_loop

        ! SUMMARIZE OUTPUT
        IF (iter > max_iter) THEN
          WRITE(*,*) 'DID NOT CONVERGE (NUMBER OF ITERATIONS:', iter, ')'
        ELSE
          WRITE(*,*) 'CONVERGED (NUMBER OF ITERATIONS:', iter, ')'
        END IF
    END SUBROUTINE solve
END MODULE subroutines


