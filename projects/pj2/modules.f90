! MAE 267
! PROJECT 1
! LOGAN HALSTROM
! 12 OCTOBER 2015

! DESCRIPTION:  Modules used for solving heat conduction of steel plate.
! Initialize and store constants used in all subroutines.

! CONTENTS:
! CONSTANTS --> Initializes constants for simulation.  Sets grid size.
! CLOCK --> Calculates clock wall-time of a process.
! MAKEGRID --> Initialize grid with correct number of points and rotation,
!                 set boundary conditions, etc.
! CELLS -->  Initialize finite volume cells and do associated calculations
! TEMPERATURE --> Calculate and store new temperature distribution
!                     for given iteration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! CONSTANTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CONSTANTS
    ! Initialize constants for simulation.  Set grid size.
    IMPLICIT NONE
    ! CFL number, for convergence (D0 is double-precision, scientific notation)
    REAL(KIND=8), PARAMETER :: CFL = 0.95D0
    ! Material constants (steel): thermal conductivity [W/(m*K)],
                                ! density [kg/m^3],
                                ! specific heat ratio [J/(kg*K)]
    REAL(KIND=8), PARAMETER :: k = 18.8D0, rho = 8000.D0, cp = 500.D0
    ! Thermal diffusivity [m^2/s]
    REAL(KIND=8), PARAMETER :: alpha = k / (cp * rho)
    ! Pi, grid rotation angle (30 deg)
    REAL(KIND=8), PARAMETER :: pi = 3.141592654D0, rot = 30.D0*pi/180.D0
    ! CPU Wall Times
    REAL(KIND=8) :: wall_time_total, wall_time_solve, wall_time_iter(1:5)
    ! read square grid size, Total grid size, size of grid on each block (local)
    INTEGER :: nx, IMAX, JMAX, IMAXBLK, JMAXBLK
    ! Dimensions of block layout, Number of Blocks,
    INTEGER :: M, N, NBLK
    ! Block boundary condition identifiers
        ! If block face is on North,east,south,west of main grid, identify which
        ! if block face is internal (has neighbor block) assign -1
    INTEGER :: NBND = 1, SBND = 2, EBND = 3, WBND = 4, NOBND = -1

CONTAINS

    SUBROUTINE read_input()
        INTEGER :: I

        ! READ INPUTS FROM FILE
            !(So I don't have to recompile each time I change an input setting)
        WRITE(*,*) 'Reading input...'
        OPEN (UNIT = 1, FILE = 'config.in')
        DO I = 1, 3
            ! Skip header lines
            READ(1,*)
        END DO
        ! READ GRIDSIZE (4th line)
        READ(1,*) nx
        ! READ BLOCKS (6th and 8th line)
        READ(1,*)
        READ(1,*) M
        READ(1,*)
        READ(1,*) N

        ! SET GRID SIZE
        IMAX = nx
        JMAX = nx
        ! CALC NUMBER OF BLOCKS
        NBLK = M * N
        ! SET SIZE OF EACH BLOCK (LOCAL MAXIMUM I, J)
        IMAXBLK = 1 + (IMAX - 1) / N
        JMAXBLK = 1 + (JMAX - 1) / M

        ! OUTPUT TO SCREEN
        WRITE(*,*) ''
        WRITE(*,*) 'Solving Mesh of size ixj:', IMAX, 'x', JMAX
        WRITE(*,*) 'With MxN blocks:', M, 'x', N
        WRITE(*,*) 'Number of blocks:', NBLK
        WRITE(*,*) ''
    END SUBROUTINE read_input
END MODULE CONSTANTS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! INITIALIZE GRID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BLOCKMOD
    ! Initialize grid with correct number of points and rotation,
    ! set boundary conditions, etc.
    USE CONSTANTS

    IMPLICIT NONE
    PUBLIC

    ! DERIVED DATA TYPE FOR GRID INFORMATION

    TYPE MESHTYPE
        ! Grid points, see cooridinate rotaion equations in problem statement
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: xp, yp, x, y
        ! Temperature at each point, temporary variable to hold temperature sum
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: T, Ttmp
        ! Iteration Parameters: timestep, cell volume, secondary cell volume,
                                    ! equation constant term
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: dt, V, V2nd, term
        ! Areas used in alternative scheme to get fluxes for second-derivative
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: Ayi, Axi, Ayj, Axj
        ! Second-derivative weighting factors for alternative distribution scheme
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: yPP, yNP, yNN, yPN
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: xNN, xPN, xPP, xNP
    END TYPE MESHTYPE

    ! DATA TYPE FOR INFORMATION ABOUT NEIGHBORS

    TYPE NBRTYPE
        INTEGER :: NBR, BC
    END TYPE NBRTYPE

    ! DERIVED DATA TYPE WITH INFORMATION PERTAINING TO SPECIFIC BLOCK

    TYPE BLKTYPE
        ! DER. DATA TYPE STORES LOCAL MESH INFO
        TYPE(MESHTYPE) :: mesh
        ! Information about face neighbors (north, east, south, west)
        TYPE(NBRTYPE) :: FaceN, FaceE, FaceS, FaceW
        ! BLOCK NUMBER
        INTEGER :: ID
        ! GLOBAL INDICIES OF MINIMUM AND MAXIMUM OF LOCAL INDICIES FOR BLOCK
        INTEGER :: IMIN, IMAX, JMIN, JMAX
    END TYPE BLKTYPE

CONTAINS
    SUBROUTINE init_blocks(blks)
        TYPE(BLKTYPE) :: blks(:)
        ! COUNTER VARIABLES
            ! IM, IN COUNT BLOCK INDICIES
            ! (IBLK COUNTS BLOCK NUMBERS, INBR IS BLOCK NEIGHBOR INDEX)
        INTEGER :: IM, IN, IBLK, INBR

        ! STEP THROUGH BLOCKS, ASSIGN IDENTIFYING INFO


        ! START AT BLOCK 1 (INCREMENT IN LOOP)
        IBLK = 0

        DO IM = 1, M
            DO IN = 1, N
                ! INCREMENT BLOCK NUMBER
                IBLK = IBLK + 1
                ! ASSIGN BLOCK NUMBER
                blks(IBLK)%ID
                ! ASSIGN GLOBAL MIN/MAX INDICIES OF LOCAL GRID
                blks(IBLK)%IMAX = 1 + IN * (IMAXBLK - 1)
                blks(IBLK)%JMAX = 1 + IM * (JMAXBLK - 1)
                blks(IBLK)%IMIN = blks(IBLK)%IMAX - (IMAXBLK - 1)
                blks(IBLK)%JMIN = blks(IBLK)%JMAX - (JMAXBLK - 1)

                ! ASSIGN NUMBERS OF FACE NEIGHBOR BLOCKS
                    !if boundary face, assign bc later
                blks(IBLK)%FaceN%NBR = IBLK + N
                blks(IBLK)%FaceS%NBR = IBLK - N
                blks(IBLK)%FaceE%NBR = IBLK + 1
                blks(IBLK)%FaceW%NBR = IBLK - 1
                ! ASSIGN FACE BOUNDARY CONDITIONS
                ! First, assign all as internal
                blks(IBLK)%FaceN%BC = -1
                blks(IBLK)%FaceS%BC = -1
                blks(IBLK)%FaceE%BC = -1
                blks(IBLK)%FaceW%BC = -1
                ! Assign faces on boundary of the actual computational grid
                ! with number corresponding to which boundary they are on
                IF ( blks(IBLK)%JMAX == JMAX ) THEN
                    ! NORTH BLOCK FACE IS ON MESH NORTH BOUNDARY
                    blks(IBLK)%FaceN%BC = NBND
                    ! un-assign neighbor that didnt really exist
                    blks(IBLK)%FaceN%NBR = -1
                ELSE IF ( blks(IBLK)%IMAX == IMAX ) THEN
                    ! EAST BLOCK FACE IS ON MESH EAST BOUNDARY
                    blks(IBLK)%FaceE%BC = EBND
                    blks(IBLK)%FaceE%NBR = -1
                ELSE IF ( blks(IBLK)%JMIN == 1 ) THEN
                    ! SOUTH BLOCK FACE IS ON MESH SOUTH BOUNDARY
                    blks(IBLK)%FaceS%BC = SBND
                    blks(IBLK)%FaceS%NBR = -1
                ELSE IF ( blks(IBLK)%IMIN == 1 ) THEN
                    ! WEST BLOCK FACE IS ON MESH WEST BOUNDARY
                    blks(IBLK)%FaceW%BC = WBND
                    blks(IBLK)%FaceW%NBR = -1
                END IF

            END DO
        END DO

    SUBROUTINE write_blocks(blks)
        ! WRITE BLOCK CONNECTIVITY FILE
        INTEGER I, BLKFILE = 99

        OPEN (UNIT = BLKFILE , FILE = "blocks.dat", form='formatted')
        DO I = 1, NBLK
            ! WRITE BLOCK NUMBER
            WRITE(BLKFILE, *), blks(I)%ID
            !
        END DO



    SUBROUTINE init_mesh(mesh)
        ! Mesh points (derived data type)
        TYPE(MESHTYPE) :: mesh
        INTEGER :: i, j

        ! ALLOCATE MESH INFORMATION
        ALLOCATE( mesh%xp(  1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%yp(  1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%x(   1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%y(   1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%T(   1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%Ttmp(1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%dt(  1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%V2nd(1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%term(1:IMAX, 1:JMAX) )

        ALLOCATE( mesh%Ayi( 1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%Axi( 1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%Ayj( 1:IMAX, 1:JMAX) )
        ALLOCATE( mesh%Axj( 1:IMAX, 1:JMAX) )

        ! ALLOCATE CELL INFORMATION
        ALLOCATE( mesh%V(  1:IMAX-1, 1:JMAX-1) )
        ALLOCATE( mesh%yPP(1:IMAX-1, 1:JMAX-1) )
        ALLOCATE( mesh%yNP(1:IMAX-1, 1:JMAX-1) )
        ALLOCATE( mesh%yNN(1:IMAX-1, 1:JMAX-1) )
        ALLOCATE( mesh%yPN(1:IMAX-1, 1:JMAX-1) )
        ALLOCATE( mesh%xNN(1:IMAX-1, 1:JMAX-1) )
        ALLOCATE( mesh%xPN(1:IMAX-1, 1:JMAX-1) )
        ALLOCATE( mesh%xPP(1:IMAX-1, 1:JMAX-1) )
        ALLOCATE( mesh%xNP(1:IMAX-1, 1:JMAX-1) )

        DO j = 1, JMAX
            DO i = 1, IMAX
                ! MAKE SQUARE GRID
                mesh%xp(i, j) = COS( 0.5D0 * pi * DFLOAT(IMAX - i) / DFLOAT(IMAX - 1) )
                mesh%yp(i, j) = COS( 0.5D0 * pi * DFLOAT(JMAX - j) / DFLOAT(JMAX - 1) )
                ! ROTATE GRID
                mesh%x(i, j) = mesh%xp(i, j) * COS(rot) + (1.D0 - mesh%yp(i, j) ) * SIN(rot)
                mesh%y(i, j) = mesh%yp(i, j) * COS(rot) + (mesh%xp(i, j)) * SIN(rot)
            END DO
        END DO
        DO j = 1, JMAX-1
            DO i = 1, IMAX-1
                ! CALC CELL VOLUME
                    ! cross product of cell diagonals p, q
                    ! where p has x,y components px, py and q likewise.
                    ! Thus, p cross q = px*qy - qx*py
                    ! where, px = x(i+1,j+1) - x(i,j), py = y(i+1,j+1) - y(i,j)
                    ! and    qx = x(i,j+1) - x(i+1,j), qy = y(i,j+1) - y(i+1,j)
                mesh%V(i,j) = ( mesh%x(i+1,j+1) - mesh%x(i,  j) ) &
                            * ( mesh%y(i,  j+1) - mesh%y(i+1,j) ) &
                            - ( mesh%x(i,  j+1) - mesh%x(i+1,j) ) &
                            * ( mesh%y(i+1,j+1) - mesh%y(i,  j) )
            END DO
        END DO
    END SUBROUTINE init_mesh

    SUBROUTINE init_temp(mesh)
        ! Initialize temperature across mesh
        ! mesh --> mesh data type
        TYPE(MESHTYPE), INTENT(INOUT) :: mesh
        INTEGER :: i, j

        !PUT DEBUG BC HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! INITIALIZE TEMPERATURE WITH DIRICHLET B.C.
        DO j = 1, JMAX
            mesh%T(1,j) = 3.D0 * mesh%yp(1,j) + 2.D0
            mesh%T(IMAX,j) = 3.D0 * mesh%yp(IMAX,j) + 2.D0
        END DO
        DO i = 1, IMAX
            mesh%T(i,1) = ABS(COS(pi * mesh%xp(i,1))) + 1.D0
            mesh%T(i,JMAX) = 5.D0 * (SIN(pi * mesh%xp(i,JMAX)) + 1.D0)
        END DO
    END SUBROUTINE init_temp

    SUBROUTINE calc_2nd_areas(m)
        ! calculate areas for secondary fluxes.
        ! m --> mesh points
        TYPE(MESHTYPE) :: m
        INTEGER :: i, j
        ! Areas used in counter-clockwise trapezoidal integration to get
        ! x and y first-derivatives for center of each cell (Green's thm)
        REAL(KIND=8) :: Ayi_half, Axi_half, Ayj_half, Axj_half

        DO j = 1, JMAX
            DO i = 1, IMAX-1
                ! CALC CELL AREAS
                m%Axj(i,j) = m%x(i+1,j) - m%x(i,j)
                m%Ayj(i,j) = m%y(i+1,j) - m%y(i,j)
            END DO
        END DO
        DO j = 1, JMAX-1
            DO i = 1, IMAX
                ! CALC CELL AREAS
                m%Axi(i,j) = m%x(i,j+1) - m%x(i,j)
                m%Ayi(i,j) = m%y(i,j+1) - m%y(i,j)
            END DO
        END DO

        ! Actual finite-volume scheme equation parameters
        DO j = 1, JMAX-1
            DO i = 1, IMAX-1

                Axi_half = ( m%Axi(i+1,j) + m%Axi(i,j) ) * 0.25D0
                Axj_half = ( m%Axj(i,j+1) + m%Axj(i,j) ) * 0.25D0
                Ayi_half = ( m%Ayi(i+1,j) + m%Ayi(i,j) ) * 0.25D0
                Ayj_half = ( m%Ayj(i,j+1) + m%Ayj(i,j) ) * 0.25D0

                ! (NN = 'negative-negative', PN = 'positive-negative',
                    ! see how fluxes are summed)
                m%xNN(i, j) = ( -Axi_half - Axj_half )
                m%xPN(i, j) = (  Axi_half - Axj_half )
                m%xPP(i, j) = (  Axi_half + Axj_half )
                m%xNP(i, j) = ( -Axi_half + Axj_half )

                m%yPP(i, j) = (  Ayi_half + Ayj_half )
                m%yNP(i, j) = ( -Ayi_half + Ayj_half )
                m%yNN(i, j) = ( -Ayi_half - Ayj_half )
                m%yPN(i, j) = (  Ayi_half - Ayj_half )
            END DO
        END DO
    END SUBROUTINE calc_2nd_areas

    SUBROUTINE calc_constants(mesh)
        ! Calculate constants for a given iteration loop.  This way,
        ! they don't need to be calculated within the loop at each iteration
        TYPE(MESHTYPE), TARGET :: mesh
        INTEGER :: i, j
        DO j = 2, JMAX - 1
            DO i = 2, IMAX - 1
                ! CALC TIMESTEP FROM CFL
                mesh%dt(i,j) = ((CFL * 0.5D0) / alpha) * mesh%V(i,j) ** 2 &
                                / ( (mesh%xp(i+1,j) - mesh%xp(i,j))**2 &
                                    + (mesh%yp(i,j+1) - mesh%yp(i,j))**2 )
                ! CALC SECONDARY VOLUMES
                ! (for rectangular mesh, just average volumes of the 4 cells
                !  surrounding the point)
                mesh%V2nd(i,j) = ( mesh%V(i,j) &
                                    + mesh%V(i-1,j) + mesh%V(i,j-1) &
                                    + mesh%V(i-1,j-1) ) * 0.25D0
                ! CALC CONSTANT TERM
                ! (this term remains constant in the equation regardless of
                !  iteration number, so only calculate once here,
                !  instead of in loop)
                mesh%term(i,j) = mesh%dt(i,j) * alpha / mesh%V2nd(i,j)
            END DO
        END DO
    END SUBROUTINE calc_constants

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! CALCULATE TEMPERATURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE calc_temp(m)
        ! Calculate first and second derivatives for finite-volume scheme
        TYPE(MESHTYPE), INTENT(INOUT) :: m
        ! First partial derivatives of temperature in x and y directions
        REAL(KIND=8) :: dTdx, dTdy
        INTEGER :: i, j

        ! RESET SUMMATION
        m%Ttmp = 0.D0

        DO j = 1, JMAX - 1
            DO i = 1, IMAX - 1
                ! CALC FIRST DERIVATIVES
                dTdx = + 0.5d0 &
                            * (( m%T(i+1,j) + m%T(i+1,j+1) ) * m%Ayi(i+1,j) &
                            -  ( m%T(i,  j) + m%T(i,  j+1) ) * m%Ayi(i,  j) &
                            -  ( m%T(i,j+1) + m%T(i+1,j+1) ) * m%Ayj(i,j+1) &
                            +  ( m%T(i,  j) + m%T(i+1,  j) ) * m%Ayj(i,  j) &
                                ) / m%V(i,j)
                dTdy = - 0.5d0 &
                            * (( m%T(i+1,j) + m%T(i+1,j+1) ) * m%Axi(i+1,j) &
                            -  ( m%T(i,  j) + m%T(i,  j+1) ) * m%Axi(i,  j) &
                            -  ( m%T(i,j+1) + m%T(i+1,j+1) ) * m%Axj(i,j+1) &
                            +  ( m%T(i,  j) + m%T(i+1,  j) ) * m%Axj(i,  j) &
                                ) / m%V(i,j)

                ! Alternate distributive scheme second-derivative operator.
                m%Ttmp(i+1,  j) = m%Ttmp(i+1,  j) + m%term(i+1,  j) * ( m%yNN(i,j) * dTdx + m%xPP(i,j) * dTdy )
                m%Ttmp(i,    j) = m%Ttmp(i,    j) + m%term(i,    j) * ( m%yPN(i,j) * dTdx + m%xNP(i,j) * dTdy )
                m%Ttmp(i,  j+1) = m%Ttmp(i,  j+1) + m%term(i,  j+1) * ( m%yPP(i,j) * dTdx + m%xNN(i,j) * dTdy )
                m%Ttmp(i+1,j+1) = m%Ttmp(i+1,j+1) + m%term(i+1,j+1) * ( m%yNP(i,j) * dTdx + m%xPN(i,j) * dTdy )
            END DO
        END DO
    END SUBROUTINE calc_temp
END MODULE BLOCKMOD


