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

MODULE CELLS
!   Initialize finite volume cells and do associated calculations
    USE MAKEGRID
    IMPLICIT NONE
    TYPE CELL
!       Cell volume
        REAL(KIND=8) :: V
!       Second-derivative weighting factors for alternative distribution scheme
        REAL(KIND=8) :: yPP, yNP, yNN, yPN
        REAL(KIND=8) :: xNN, xPN, xPP, xNP
    END TYPE CELL

CONTAINS
    SUBROUTINE init_cells(cells, mesh)
!       cells --> derived data type containing cell info
!       mesh --> derived data type containing mesh point info
        TYPE(CELL), TARGET :: cells(1:IMAX-1,1:JMAX-1)
        TYPE(GRID) :: mesh(1:IMAX, 1:JMAX)
        INTEGER :: i, j
        DO J = 1, JMAX
            DO I = 1, IMAX
!               CALC CELL VOLUMES
!                   (length in x-dir times length in y-dir)
                cells(i,j)%V = (mesh(i+1,j)%xp - mesh(i,j)%xp) &
                                    * (mesh(i,j+1).yp - mesh(i,j).yp)
            END DO
        END DO
    END SUBROUTINE init_cells

    SUBROUTINE calc_2nd_areas(cells, m)
!       calculate areas for secondary fluxes.
!       cells --> derived data type with cell data, target for c
!       m --> mesh points
        TYPE(GRID), TARGET :: m(1:IMAX, 1:JMAX)
        TYPE(CELL), TARGET :: cells(1:IMAX-1, 1:JMAX-1)
        TYPE(CELL), POINTER :: c
        INTEGER :: i, j
!       Areas used in alternative scheme to get fluxes for second-derivative
        REAL(KIND=8) :: Ayi, Axi, Ayj, Axj
!       Areas used in counter-clockwise trapezoidal integration to get
!       x and y first-derivatives for center of each cell (Green's thm)
        REAL(KIND=8) :: Ayi_half, Axi_half, Ayj_half, Axj_half

!       CALC CELL AREAS
        Axi(i,j) = m(i,j+1)%x - m(i,j)%x
        Axj(i,j) = m(i+1,j)%x - m(i,j)%x
        Ayi(i,j) = m(i,j+1)%y - m(i,j)%y
        Ayj(i,j) = m(i+1,j)%y - m(i,j)%y

        Axi_half(i,j) = ( Axi(i+1,j) + Axi(i,j) ) * 0.25D0
        Axj_half(i,j) = ( Axj(i,j+1) + Axj(i,j) ) * 0.25D0
        Ayi_half(i,j) = ( Ayi(i+1,j) + Ayi(i,j) ) * 0.25D0
        Ayj_half(i,j) = ( Ayj(i,j+1) + Ayj(i,j) ) * 0.25D0

!       Actual finite-volume scheme equation parameters
        DO J = 1, JMAX
            DO I = 1, IMAX
!               (NN = 'negative-negative', PN = 'positive-negative',
!                   see how fluxes are summed)
                c%xNN = ( -Axi_half(i,j) - Axj_half(i,j) )
                c%xPN = (  Axi_half(i,j) - Axj_half(i,j) )
                c%xPP = (  Axi_half(i,j) + Axj_half(i,j) )
                c%xNP = ( -Axi_half(i,j) + Axj_half(i,j) )

                c%yPP = (  Ayi_half(i,j) + Ayj_half(i,j) )
                c%yNP = ( -Ayi_half(i,j) + Ayj_half(i,j) )
                c%yNN = ( -Ayi_half(i,j) - Ayj_half(i,j) )
                c%yPN = (  Ayi_half(i,j) - Ayj_half(i,j) )
            END DO
        END DO
    END SUBROUTINE calc_2nd_areas


