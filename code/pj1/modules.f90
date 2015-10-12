! MAE 267
! PROJECT 1
! LOGAN HALSTROM
! 12 OCTOBER 2015

! DESCRIPTION:  Modules used for solving heat conduction of steel plate.
! Initialize and store constants used in all subroutines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! CONSTANTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CONSTANTS
!   Initialize constants for simulation.  Set grid size.
    IMPLICIT NONE
!   CFL number, for convergence (D0 is double-precision, scientific notation)
    REAL(KIND=8), PARAMETER :: CFL = 0.5D0
!   Material constants (steel): thermal conductivity [W/(m*K)],
!                               density [kg/m^3],
!                               specific heat ratio [J/(kg*K)]
    REAL(KIND=8), PARAMETER :: k = 18.8D0, rho = 8000.D0, cp = 500.D0
!   Thermal diffusivity [m^2/s]
    REAL(KIND=8), PARAMETER :: alpha = k / (cp * rho)
!   Pi, grid rotation angle (30 deg)
    REAL(KIND=8), PARAMETER :: pi = 3.141592654D0, rot = 30.D0*pi/180.D0
!   Grid size
    INTEGER :: IMAX, JMAX

CONTAINS
    SUBROUTINE GRIDSIZE(n)
!           Set size of grid (square)
        INTEGER :: n
        IMAX = N
        JMAX = N
    END SUBROUTINE GRIDSIZE

END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! WALL CLOCK TIME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CLOCK
!   Clock wall-time of a process.
    INTEGER clock_start, clock_end, clock_max, clock_rate
    REAL(KIND=8) wall_time

CONTAINS
    SUBROUTINE start_clock()
!       get clock parameters
        CALL SYSTEM_CLOCK(count_max=clock_max, count_rate=clock_rate)
!       Get start time
        CALL SYSTEM_CLOCK(clock_start)
    END SUBROUTINE start_clock()

    SUBROUTINE end_clock()
!       Get end time
        CALL SYSTEM_CLOCK(clock_end)
        wall_time = DFLOAT(clock_end - clock_start) / DFLOAT(clock_rate)
        PRINT*, 'Solver wall clock time (seconds):', wall_time
    END SUBROUTINE end_clock
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! INITIALIZE GRID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MAKEGRID
!   Initialize grid with correct number of points and rotation,
!   set boundary conditions, etc.
    USE CONSTANTS
    IMPLICIT NONE

    PUBLIC
!   Derived data type
    TYPE GRID
        INTEGER :: i, j
!       #Grid points, see cooridinate rotaion equations in problem statement
        REAL(KIND=8) :: xp, yp, x, y
!       #Temperature at each point
        REAL(KIND=8) :: T

!        real(kind=8) :: timestep, Vol2, const!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END TYPE GRID

CONTAINS
    SUBROUTINE init_grid(mesh)
!       Mesh points (derived data type)
        TYPE(GRID), TARGET :: mesh(1:IMAX, 1:JMAX)
!       Pointer for mesh points
        TYPE(GRID), POINTER :: p
        INTEGER :: i, j

        DO J = 1, JMAX
            DO I = 1, IMAX
                p => mesh(i, j)
!               'p' points to 'mesh', i is variable in derived data type
!                   accessed by '%'
!               MAKE SQUARE GRID
                p%i = i
                p%j = j
!               ROTATE GRID
                p%xp = COS( 0.5D0 * pi * DFLOAT(IMAX - i) / DFLOAT(IMAX - 1) )
                p%yp = COS( 0.5D0 * pi * DFLOAT(JMAX - j) / DFLOAT(JMAX - 1) )

                p%x = p%xp * COS(rot) + (1.D0 - p%yp ) * SIN(rot)
                p%y = p%yp * COS(rot) + (p%xp) * SIN(rot)
            END DO
        END DO
    END SUBROUTINE init_grid

    SUBROUTINE init_temp(p, T)
!       Initialize temperature across mesh
!       p --> pointer for mesh vector
!       T --> initial temperature profile
        TYPE(GRID), INTENT(INOUT) :: p
        REAL(KIND=8) :: T
!       SET MESH POINTS WITH INITIAL TEMPERATURE PROFILE
        p%T = T
    END SUBROUTINE init_temp
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! CELLS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




