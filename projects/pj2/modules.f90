! MAE 267
! PROJECT 2
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
                                ! initial temperature
    REAL(KIND=8), PARAMETER :: k = 18.8D0, rho = 8000.D0, cp = 500.D0, T0 = 3.5D0
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
        ! If block face is on North,east,south,west of main grid, identify
    INTEGER :: NBND = 1, SBND = 2, EBND = 3, WBND = 4
    ! Output directories
!     CHARACTER(LEN=7) :: outdir = "Results", casedir
    CHARACTER(LEN=16) :: casedir

CONTAINS

    SUBROUTINE read_input()
        INTEGER :: I
        CHARACTER(LEN=3) :: strNX
        CHARACTER(LEN=1) :: strN, strM
        LOGICAL :: exists

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

        ! OUTPUT DIRECTORIES
        ! write integers to strings
        WRITE(strNX, '(I3)') nx
        WRITE(strN,  '(I1)') N
        WRITE(strM,  '(I1)') M
        !case output directory: nx_NxM (i.e. 'Results/101_5x4')
        casedir = 'Results/' // strNX // '_' // strN // 'x' // strM // '/'
        ! MAKE DIRECTORIES (IF THEY DONT ALREADY EXIST)
!         CALL EXECUTE_COMMAND_LINE ("mkdir -p " // outdir // '/' // casedir)
        CALL EXECUTE_COMMAND_LINE ("mkdir -p " // casedir)

        ! OUTPUT TO SCREEN
        WRITE(*,*) ''
        WRITE(*,*) 'Solving Mesh of size ixj:', IMAX, 'x', JMAX
        WRITE(*,*) 'With MxN blocks:', M, 'x', N
        WRITE(*,*) 'Number of blocks:', NBLK
        WRITE(*,*) 'Block size ixj:', IMAXBLK, 'x', JMAXBLK
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
        ! STORE NUMBER OF NEIGHBOR BLOCK, OR NUMBER OF BOUNDARY CONDITION
        INTEGER :: BC, NB
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
        ! BLOCK ORIENTATION
        INTEGER :: orientation
    END TYPE BLKTYPE

CONTAINS
    SUBROUTINE init_blocks(b)
        ! BLOCK DATA TYPE
        TYPE(BLKTYPE) :: b(:)
        ! COUNTER VARIABLES
            ! IM, IN COUNT BLOCK INDICIES
            ! (IBLK COUNTS BLOCK NUMBERS, INBR IS BLOCK NEIGHBOR INDEX)
        INTEGER :: I, J, IBLK, INBR

        ! STEP THROUGH BLOCKS, ASSIGN IDENTIFYING INFO


        ! START AT BLOCK 1 (INCREMENT IN LOOP)
        IBLK = 0

        DO J = 1, M
            DO I = 1, N
                ! INCREMENT BLOCK NUMBER
                IBLK = IBLK + 1
                ! ASSIGN BLOCK NUMBER
                b(IBLK)%ID = IBLK
                ! ASSIGN GLOBAL MIN/MAX INDICIES OF LOCAL GRID
                b(IBLK)%IMIN = 1 + (IMAXBLK - 1) * (I - 1)
                b(IBLK)%JMIN = 1 + (JMAXBLK - 1) * (J - 1)
                b(IBLK)%IMAX = b(IBLK)%IMIN + (IMAXBLK - 1)
                b(IBLK)%JMAX = b(IBLK)%JMIN + (JMAXBLK - 1)

                ! ASSIGN NUMBERS OF FACE NEIGHBOR BLOCKS
                    !if boundary face, assign bc later
                b(IBLK)%FaceN%NB = IBLK + N
                b(IBLK)%FaceS%NB = IBLK - N
                b(IBLK)%FaceE%NB = IBLK + 1
                b(IBLK)%FaceW%NB = IBLK - 1

                ! ASSIGN FACE BOUNDARY CONDITIONS
                    ! initialize as all internal
                b(IBLK)%FaceN%BC = -1
                b(IBLK)%FaceS%BC = -1
                b(IBLK)%FaceE%BC = -1
                b(IBLK)%FaceW%BC = -1
                ! Assign faces on boundary of the actual computational grid
                ! with number corresponding to which boundary they are on
                IF ( b(IBLK)%JMAX == JMAX ) THEN
                    ! NORTH BLOCK FACE IS ON MESH NORTH BOUNDARY
                    b(IBLK)%FaceN%BC = NBND
                    ! un-assign neighbor that wasnt real
                    b(IBLK)%FaceN%NB = -1
                END IF
                IF ( b(IBLK)%IMAX == IMAX ) THEN
                    ! EAST BLOCK FACE IS ON MESH EAST BOUNDARY
                    b(IBLK)%FaceE%BC = EBND
                    b(IBLK)%FaceE%NB = -1
                END IF
                IF ( b(IBLK)%JMIN == 1 ) THEN
                    ! SOUTH BLOCK FACE IS ON MESH SOUTH BOUNDARY
                    b(IBLK)%FaceS%BC = SBND
                    b(IBLK)%FaceS%NB = -1
                END IF
                IF ( b(IBLK)%IMIN == 1 ) THEN
                    ! WEST BLOCK FACE IS ON MESH WEST BOUNDARY
                    b(IBLK)%FaceW%BC = WBND
                    b(IBLK)%FaceW%NB = -1
                END IF

                ! BLOCK ORIENTATION
                    ! same for all in this project
                b(IBLK)%orientation = 1

            END DO
        END DO
    END SUBROUTINE init_blocks

    SUBROUTINE write_blocks(b)
        ! WRITE BLOCK CONNECTIVITY FILE

        ! BLOCK DATA TYPE
        TYPE(BLKTYPE) :: b(:)
        INTEGER :: I, BLKFILE = 99

        11 format(3I5)
        22 format(33I5)

        OPEN (UNIT = BLKFILE , FILE = casedir // "blockconfig.dat", form='formatted')
        ! WRITE AMOUNT OF BLOCKS AND DIMENSIONS
        WRITE(BLKFILE, 11) NBLK, IMAXBLK, JMAXBLK
        DO I = 1, NBLK
            ! FOR EACH BLOCK, WRITE BLOCK NUMBER, STARTING/ENDING GLOBAL INDICES.
            ! THEN BOUNDARY CONDITION AND NEIGHBOR NUMBER FOR EACH FACE:
            ! NORTH EAST SOUTH WEST
            WRITE(BLKFILE, 22) b(I)%ID, &
                b(I)%IMIN, b(I)%JMIN, &
                b(I)%IMAX, b(I)%JMAX, &
                b(I)%FaceN%BC, b(I)%FaceN%NB, &
                b(I)%FaceE%BC, b(I)%FaceE%NB, &
                b(I)%FaceS%BC, b(I)%FaceS%NB, &
                b(I)%FaceW%BC, b(I)%FaceW%NB, &
                b(I)%orientation
        END DO
        CLOSE(BLKFILE)
    END SUBROUTINE write_blocks

    SUBROUTINE init_mesh(b)
        ! BLOCK DATA TYPE
        TYPE(BLKTYPE) :: b(:)
        INTEGER :: IBLK, I, J

        DO IBLK = 1, NBLK

            ! ALLOCATE MESH INFORMATION
                ! ADD EXTRA INDEX AT BEGINNING AND END FOR GHOST NODES
            ALLOCATE( b(IBLK)%mesh%xp(  0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%yp(  0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%x(   0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%y(   0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%T(   0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%Ttmp(0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%dt(  0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%V2nd(0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%term(0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%Ayi( 0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%Axi( 0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%Ayj( 0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%Axj( 0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( b(IBLK)%mesh%V(   0:IMAXBLK-1+1, 0:JMAXBLK-1+1) )
            ALLOCATE( b(IBLK)%mesh%yPP( 0:IMAXBLK-1+1, 0:JMAXBLK-1+1) )
            ALLOCATE( b(IBLK)%mesh%yNP( 0:IMAXBLK-1+1, 0:JMAXBLK-1+1) )
            ALLOCATE( b(IBLK)%mesh%yNN( 0:IMAXBLK-1+1, 0:JMAXBLK-1+1) )
            ALLOCATE( b(IBLK)%mesh%yPN( 0:IMAXBLK-1+1, 0:JMAXBLK-1+1) )
            ALLOCATE( b(IBLK)%mesh%xNN( 0:IMAXBLK-1+1, 0:JMAXBLK-1+1) )
            ALLOCATE( b(IBLK)%mesh%xPN( 0:IMAXBLK-1+1, 0:JMAXBLK-1+1) )
            ALLOCATE( b(IBLK)%mesh%xPP( 0:IMAXBLK-1+1, 0:JMAXBLK-1+1) )
            ALLOCATE( b(IBLK)%mesh%xNP( 0:IMAXBLK-1+1, 0:JMAXBLK-1+1) )

            ! STEP THROUGH LOCAL INDICIES OF EACH BLOCK
            DO J = 0, JMAXBLK+1
                DO I = 0, IMAXBLK+1
                    ! MAKE SQUARE GRID
                        ! CONVERT FROM LOCAL TO GLOBAL INDEX:
                            ! Iglobal = Block%IMIN + (Ilocal - 1)
                    b(IBLK)%mesh%xp(I, J) = COS( 0.5D0 * pi * DFLOAT(IMAX - ( b(IBLK)%IMIN + I - 1) ) / DFLOAT(IMAX - 1) )
                    b(IBLK)%mesh%yp(I, J) = COS( 0.5D0 * pi * DFLOAT(JMAX - ( b(IBLK)%JMIN + J - 1) ) / DFLOAT(JMAX - 1) )
                    ! ROTATE GRID
                    b(IBLK)%mesh%x(I, J) = b(IBLK)%mesh%xp(I, J) * COS(rot) + (1.D0 - b(IBLK)%mesh%yp(I, J) ) * SIN(rot)
                    b(IBLK)%mesh%y(I, J) = b(IBLK)%mesh%yp(I, J) * COS(rot) + (       b(IBLK)%mesh%xp(I, J) ) * SIN(rot)
                END DO
            END DO
            DO J = 0, JMAXBLK+1-1
                DO I = 0, IMAXBLK+1-1
                    ! CALC CELL VOLUME
                        ! cross product of cell diagonals p, q
                        ! where p has x,y components px, py and q likewise.
                        ! Thus, p cross q = px*qy - qx*py
                        ! where, px = x(i+1,j+1) - x(i,j), py = y(i+1,j+1) - y(i,j)
                        ! and    qx = x(i,j+1) - x(i+1,j), qy = y(i,j+1) - y(i+1,j)
                    b(IBLK)%mesh%V(I,J) = ( b(IBLK)%mesh%x(I+1,J+1) &
                                          - b(IBLK)%mesh%x(I,  J) ) &
                            * ( b(IBLK)%mesh%y(I,  J+1) - b(IBLK)%mesh%y(I+1,J) ) &
                            - ( b(IBLK)%mesh%x(I,  J+1) - b(IBLK)%mesh%x(I+1,J) ) &
                            * ( b(IBLK)%mesh%y(I+1,J+1) - b(IBLK)%mesh%y(I,  J) )
                END DO
            END DO
        END DO
    END SUBROUTINE init_mesh

    SUBROUTINE init_temp(blocks)
        ! Initialize temperature across mesh
        ! BLOCK DATA TYPE
        TYPE(BLKTYPE), TARGET  :: blocks(:)
        TYPE(BLKTYPE), POINTER :: b
        TYPE(MESHTYPE), POINTER :: m
        INTEGER :: IBLK, I, J


        !PUT DEBUG BC HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO IBLK = 1, NBLK
            b => blocks(IBLK)
            m => blocks(IBLK)%mesh
            ! FIRST, INITIALIZE ALL POINT TO INITIAL TEMPERATURE (T0)
            m%T(1:IMAXBLK, 1:JMAXBLK) = T0
            ! THEN, INITIALIZE BOUNDARIES DIRICHLET B.C.
            ! face on north boundary
            IF (b%FaceN%BC == NBND) THEN
                DO I = 1, IMAXBLK
                    m%T(I,JMAX) = 5.D0 * (SIN(pi * m%xp(I,JMAX)) + 1.D0)
                END DO
            END IF
            IF (b%FaceS%BC == SBND) THEN
                DO I = 1, IMAXBLK
                    m%T(I,1) = ABS(COS(pi * m%xp(I,1))) + 1.D0
                END DO
            END IF
            IF (b%FaceE%BC == EBND) THEN
                DO J = 1, JMAXBLK
                    m%T(IMAX,J) = 3.D0 * m%yp(IMAX,J) + 2.D0
                END DO
            END IF
            IF (b%FaceW%BC == WBND) THEN
                DO J = 1, JMAXBLK
                    m%T(I,1) = ABS(COS(pi * m%xp(I,1))) + 1.D0
                END DO
            END IF

!             DO J = 1, JMAXBLK
!                 b(IBLK)%mesh%T(1,J) = 3.D0 * b(IBLK)%mesh%yp(1,J) + 2.D0
!                 b(IBLK)%mesh%T(IMAX,J) = 3.D0 * b(IBLK)%mesh%yp(IMAX,J) + 2.D0
!             END DO
!             DO I = 1, IMAXBLK
!                 b(IBLK)%mesh%T(I,1) = ABS(COS(pi * b(IBLK)%mesh%xp(I,1))) + 1.D0
!                 b(IBLK)%mesh%T(I,JMAX) = 5.D0 * (SIN(pi * b(IBLK)%mesh%xp(I,JMAX)) + 1.D0)
!             END DO
        END DO
    END SUBROUTINE init_temp

    SUBROUTINE calc_2nd_areas(blocks)
        ! calculate areas for secondary fluxes.
        ! BLOCK DATA TYPE
        TYPE(BLKTYPE), TARGET :: blocks(:)
        TYPE(MESHTYPE), POINTER :: m
        INTEGER :: IBLK, I, J
        ! Areas used in counter-clockwise trapezoidal integration to get
        ! x and y first-derivatives for center of each cell (Green's thm)
        REAL(KIND=8) :: Ayi_half, Axi_half, Ayj_half, Axj_half

        DO IBLK = 1, NBLK
            m => blocks(IBLK)%mesh
            DO J = 1, JMAX
                DO I = 1, IMAX-1
                    ! CALC CELL AREAS
                    m%Axj(I,J) = m%x(I+1,J) - m%x(I,J)
                    m%Ayj(I,J) = m%y(I+1,J) - m%y(I,J)
                END DO
            END DO
            DO J = 1, JMAX-1
                DO I = 1, IMAX
                    ! CALC CELL AREAS
                    m%Axi(I,J) = m%x(I,J+1) - m%x(I,J)
                    m%Ayi(I,J) = m%y(I,J+1) - m%y(I,J)
                END DO
            END DO

            ! Actual finite-volume scheme equation parameters
            DO J = 1, JMAX-1
                DO I = 1, IMAX-1

                    Axi_half = ( m%Axi(I+1,J) + m%Axi(I,J) ) * 0.25D0
                    Axj_half = ( m%Axj(I,J+1) + m%Axj(I,J) ) * 0.25D0
                    Ayi_half = ( m%Ayi(I+1,J) + m%Ayi(I,J) ) * 0.25D0
                    Ayj_half = ( m%Ayj(I,J+1) + m%Ayj(I,J) ) * 0.25D0

                    ! (NN = 'negative-negative', PN = 'positive-negative',
                        ! see how fluxes are summed)
                    m%xNN(I, J) = ( -Axi_half - Axj_half )
                    m%xPN(I, J) = (  Axi_half - Axj_half )
                    m%xPP(I, J) = (  Axi_half + Axj_half )
                    m%xNP(I, J) = ( -Axi_half + Axj_half )
                    m%yPP(I, J) = (  Ayi_half + Ayj_half )
                    m%yNP(I, J) = ( -Ayi_half + Ayj_half )
                    m%yNN(I, J) = ( -Ayi_half - Ayj_half )
                    m%yPN(I, J) = (  Ayi_half - Ayj_half )
                END DO
            END DO
        END DO
    END SUBROUTINE calc_2nd_areas

    SUBROUTINE calc_constants(blocks)
        ! Calculate constants for a given iteration loop.  This way,
        ! they don't need to be calculated within the loop at each iteration
        TYPE(BLKTYPE), TARGET :: blocks(:)
        TYPE(MESHTYPE), POINTER :: m
        INTEGER :: IBLK, I, J
        DO IBLK = 1, NBLK
            m => blocks(IBLK)%mesh
            DO J = 2, JMAX - 1
                DO I = 2, IMAX - 1
                    ! CALC TIMESTEP FROM CFL
                    m%dt(I,J) = ((CFL * 0.5D0) / alpha) * m%V(I,J) ** 2 &
                                    / ( (m%xp(I+1,J) - m%xp(I,J))**2 &
                                        + (m%yp(I,J+1) - m%yp(I,J))**2 )
                    ! CALC SECONDARY VOLUMES
                    ! (for rectangular mesh, just average volumes of the 4 cells
                    !  surrounding the point)
                    m%V2nd(I,J) = ( m%V(I,J) &
                                        + m%V(I-1,J) + m%V(I,J-1) &
                                        + m%V(I-1,J-1) ) * 0.25D0
                    ! CALC CONSTANT TERM
                    ! (this term remains constant in the equation regardless of
                    !  iteration number, so only calculate once here,
                    !  instead of in loop)
                    m%term(I,J) = m%dt(I,J) * alpha / m%V2nd(I,J)
                END DO
            END DO
        END DO
    END SUBROUTINE calc_constants

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! CALCULATE TEMPERATURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE calc_temp(blocks)
        ! Calculate first and second derivatives for finite-volume scheme
        TYPE(BLKTYPE), TARGET :: blocks(:)
        TYPE(MESHTYPE), POINTER :: m
        ! First partial derivatives of temperature in x and y directions
        REAL(KIND=8) :: dTdx, dTdy
        INTEGER :: IBLK, I, J

        DO IBLK = 1, NBLK
            m => blocks(IBLK)%mesh

            ! RESET SUMMATION
            m%Ttmp = 0.D0

            DO J = 1, JMAX - 1
                DO I = 1, IMAX - 1
                    ! CALC FIRST DERIVATIVES
                    dTdx = + 0.5d0 &
                                * (( m%T(I+1,J) + m%T(I+1,J+1) ) * m%Ayi(I+1,J) &
                                -  ( m%T(I,  J) + m%T(I,  J+1) ) * m%Ayi(I,  J) &
                                -  ( m%T(I,J+1) + m%T(I+1,J+1) ) * m%Ayj(I,J+1) &
                                +  ( m%T(I,  J) + m%T(I+1,  J) ) * m%Ayj(I,  J) &
                                    ) / m%V(I,J)
                    dTdy = - 0.5d0 &
                                * (( m%T(I+1,J) + m%T(I+1,J+1) ) * m%Axi(I+1,J) &
                                -  ( m%T(I,  J) + m%T(I,  J+1) ) * m%Axi(I,  J) &
                                -  ( m%T(I,J+1) + m%T(I+1,J+1) ) * m%Axj(I,J+1) &
                                +  ( m%T(I,  J) + m%T(I+1,  J) ) * m%Axj(I,  J) &
                                    ) / m%V(I,J)

                    ! Alternate distributive scheme second-derivative operator.
                    m%Ttmp(I+1,  J) = m%Ttmp(I+1,  J) + m%term(I+1,  J) * ( m%yNN(I,J) * dTdx + m%xPP(I,J) * dTdy )
                    m%Ttmp(I,    J) = m%Ttmp(I,    J) + m%term(I,    J) * ( m%yPN(I,J) * dTdx + m%xNP(I,J) * dTdy )
                    m%Ttmp(I,  J+1) = m%Ttmp(I,  J+1) + m%term(I,  J+1) * ( m%yPP(I,J) * dTdx + m%xNN(I,J) * dTdy )
                    m%Ttmp(I+1,J+1) = m%Ttmp(I+1,J+1) + m%term(I+1,J+1) * ( m%yNP(I,J) * dTdx + m%xPN(I,J) * dTdy )
                END DO
            END DO
        END DO
    END SUBROUTINE calc_temp
END MODULE BLOCKMOD


