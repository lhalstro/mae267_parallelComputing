! MAE 267
! PROJECT 4
! LOGAN HALSTROM
! 14 NOVEMBER 2015

! DESCRIPTION:  This module contains functions for information input and output.
! Write grid and temperature files in PLOT3D format.
! Write and read block grid configuration file

! NOTE: How to Visualize Blocks in Paraview:
    ! open unformatted PLOT3D file.
    ! Change 'Coloring' from 'Solid' to 'vtkCompositeIndex'

MODULE IO
    USE CONSTANTS
    USE BLOCKMOD
    IMPLICIT NONE

    ! VARIABLES
    INTEGER :: gridUnit  = 30   ! Unit for grid file
    INTEGER :: tempUnit = 21    ! Unit for temp file
    INTEGER :: resUnit = 23
    REAL(KIND=8) :: tRef = 1.D0          ! tRef number
    REAL(KIND=8) :: dum = 0.D0          ! dummy values

    ! LINKED LIST OF RESIDUAL HISTORY

    TYPE RESLIST
        ! Next element in linked list
        TYPE(RESLIST), POINTER :: next
        ! items in link:
        REAL(KIND=8) :: res
        INTEGER :: iter
    END TYPE RESLIST

    CONTAINS

        SUBROUTINE write_config(procs)
        ! Write block connectivity file with neighbor and BC info
        ! for each processor.
        ! Also write PLOT3D restart files for each processor.

        TYPE(PROCTYPE), TARGET :: procs(:)
        TYPE(PROCTYPE), POINTER :: p
        ! BLOCK DATA TYPE
        TYPE(BLKTYPE), POINTER :: b
        INTEGER :: IP, IB, BLKFILE = 99
        CHARACTER(2) :: procname
        CHARACTER(20) :: xfile, qfile

        33 FORMAT(A)
        11 FORMAT( 3I7)
        22 FORMAT(33I7)
        44 FORMAT(33A7)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! WRITE CONFIG FILE FOR EACH PROCESSOR !!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO IP = 1, NPROCS
            p => procs(IP)

            ! MAKE FILE NAME (i.e. 'p01.config')
            IF (p%ID<10) THEN
                ! IF SINGLE DIGIT, PAD WITH 0 IN FRONT
                WRITE(procname, '(A,I1)') '0', p%ID
            ELSE
                WRITE(procname, '(I2)') p%ID
            END IF

            OPEN (UNIT = BLKFILE , FILE = TRIM("p" // procname // ".config"), form='formatted')

            ! WRITE AMOUNT OF BLOCKS AND DIMENSIONS
            WRITE(BLKFILE, 33) 'NBLK' // ' IMAXBLK' // ' JMAXBLK'
            WRITE(BLKFILE, 11) p%NBLK, IMAXBLK, JMAXBLK

            ! HEADER
            WRITE(BLKFILE, 44) 'ID', 'IMIN', 'JMIN', 'SIZE', &
                               'IMNL', 'IMXL', 'JMNL', 'JMXL', &
                               'NNB', 'NNP', 'NLOC', &
                               'SNB', 'SNP', 'SLOC', &
                               'ENB', 'ENP', 'ELOC', &
                               'WNB', 'WNP', 'WLOC', &
                               'NENB', 'NENP', 'NEL', &
                               'SENB', 'SENP', 'SEL', &
                               'SWNB', 'SWNP', 'SWL', &
                               'NWNB', 'NWNP', 'NWL', &
                               'ORI'
            DO IB = 1, p%NBLK
                b => p%blocks(IB)
                ! FOR EACH BLOCK, WRITE BLOCK NUMBER, STARTING/ENDING GLOBAL INDICES.
                ! THEN BOUNDARY CONDITION AND NEIGHBOR NUMBER FOR EACH FACE:
                ! NORTH EAST SOUTH WEST
                WRITE(BLKFILE, 22) b%ID, b%IMIN, b%JMIN, INT(b%SIZE), &
                                   b%IMINLOC, b%IMAXLOC, b%JMINLOC, b%JMAXLOC, &
                                   b%NB%N,  b%NP%N,  b%NBLOC%N, &
                                   b%NB%S,  b%NP%S,  b%NBLOC%S, &
                                   b%NB%E,  b%NP%E,  b%NBLOC%E, &
                                   b%NB%W,  b%NP%W,  b%NBLOC%W, &
                                   b%NB%NE, b%NP%NE, b%NBLOC%NE, &
                                   b%NB%SE, b%NP%SE, b%NBLOC%SE, &
                                   b%NB%SW, b%NP%SW, b%NBLOC%SW, &
                                   b%NB%NW, b%NP%NW, b%NBLOC%NW, &
                                   b%ORIENT
            END DO
            CLOSE(BLKFILE)
        END DO

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! WRITE SOLUTION RESTART FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO IP = 1, NPROCS
            p => procs(IP)
            ! MAKE FILE NAME
            IF (p%ID<10) THEN
                ! IF SINGLE DIGIT, PAD WITH 0 IN FRONT
                WRITE(procname, '(A,I1)') '0', p%ID
            ELSE
                WRITE(procname, '(I2)') p%ID
            END IF
            xfile = "p" // procname // ".grid"
            qfile = "p" // procname // ".T"
            CALL plot3D(p%blocks, p%NBLK, xfile, qfile)
        END DO
    END SUBROUTINE write_config

    SUBROUTINE read_config(blocks)
        ! Called by each processor individually for its own blocks
        ! For given processor, read corresponding configuration file.
        ! Get neighbor connectivity info
        ! Also read PLOT3D restart files for grids for given processor

        ! BLOCK DATA TYPE
        TYPE(BLKTYPE), POINTER :: blocks(:)
        TYPE(BLKTYPE), POINTER :: b
        TYPE(MESHTYPE), POINTER :: m
        INTEGER :: IP, IB, BLKFILE = 99
        CHARACTER(2) :: procname
        CHARACTER(20) :: xfile, qfile

        33 FORMAT(A)
        11 FORMAT( 3I7)
        22 FORMAT(33I7)
        44 FORMAT(33A7)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! READ CONFIG FILE FOR GIVEN PRO!SSOR !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! FILE NAME (i.e. 'p01.config')
        IF (MYID<10) THEN
            ! IF SINGLE DIGIT, PAD WITH 0 IN FRONT
            WRITE(procname, '(A,I1)') '0', MYID
        ELSE
            WRITE(procname, '(I2)') MYID
        END IF

        OPEN (UNIT = BLKFILE , FILE = TRIM("p" // procname // ".config"), form='formatted')

        ! WRITE AMOUNT OF BLOCKS AND DIMENSIONS
        READ(BLKFILE, *)
        READ(BLKFILE, 11) MYNBLK, IMAXBLK, JMAXBLK
        ! (MYNBLK is global variable on this processor that is the number of
        !   blocks allocated to this processor)

        ! ALLOCATE BLOCKS FOR THIS PROCESSOR
        ! ('blocks' stores just the blocks for this processor.  It is
        !  in parallel for each processor)
        ALLOCATE( blocks(1:MYNBLK) )

        ! ALLOCATE MESH STUFF TO BE READ IN
        DO IB = 1, MYNBLK
            m => blocks(IB)%mesh

            ! ALLOCATE MESH INFORMATION
                ! ADD EXTRA INDEX AT BEGINNING AND END FOR GHOST NODES
            ALLOCATE( m%xp(  0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%yp(  0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%x(   0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%y(   0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%T(   0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%Ttmp(0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%dt(  0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%V2nd(0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%term(0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%Ayi( 0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%Axi( 0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%Ayj( 0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%Axj( 0:IMAXBLK+1,   0:JMAXBLK+1) )
            ALLOCATE( m%V(   0:IMAXBLK,     0:JMAXBLK  ) )
            ALLOCATE( m%yPP( 0:IMAXBLK,     0:JMAXBLK  ) )
            ALLOCATE( m%yNP( 0:IMAXBLK,     0:JMAXBLK  ) )
            ALLOCATE( m%yNN( 0:IMAXBLK,     0:JMAXBLK  ) )
            ALLOCATE( m%yPN( 0:IMAXBLK,     0:JMAXBLK  ) )
            ALLOCATE( m%xNN( 0:IMAXBLK,     0:JMAXBLK  ) )
            ALLOCATE( m%xPN( 0:IMAXBLK,     0:JMAXBLK  ) )
            ALLOCATE( m%xPP( 0:IMAXBLK,     0:JMAXBLK  ) )
            ALLOCATE( m%xNP( 0:IMAXBLK,     0:JMAXBLK  ) )
        END DO

        ! HEADER
        READ(BLKFILE, *)
        DO IB = 1, MYNBLK
            b => blocks(IB)
            ! FOR EACH BLOCK, READ BLOCK NUMBER, STARTING/ENDING GLOBAL INDICES.
            ! THEN BOUNDARY CONDITION AND NEIGHBOR NUMBER FOR EACH FACE:
            ! NORTH EAST SOUTH WEST
            READ(BLKFILE, 22) b%ID, b%IMIN, b%JMIN, b%SIZE, &
                              b%IMINLOC, b%IMAXLOC, b%JMINLOC, b%JMAXLOC, &
                              b%NB%N,  b%NP%N,  b%NBLOC%N, &
                              b%NB%S,  b%NP%S,  b%NBLOC%S, &
                              b%NB%E,  b%NP%E,  b%NBLOC%E, &
                              b%NB%W,  b%NP%W,  b%NBLOC%W, &
                              b%NB%NE, b%NP%NE, b%NBLOC%NE, &
                              b%NB%SE, b%NP%SE, b%NBLOC%SE, &
                              b%NB%SW, b%NP%SW, b%NBLOC%SW, &
                              b%NB%NW, b%NP%NW, b%NBLOC%NW, &
                              b%ORIENT
        END DO
        CLOSE(BLKFILE)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! READ SOLUTION RESTART FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! MAKE FILE NAME
        IF (MYID<10) THEN
            ! IF SINGLE DIGIT, PAD WITH 0 IN FRONT
            WRITE(procname, '(A,I1)') '0', MYID
        ELSE
            WRITE(procname, '(I2)') MYID
        END IF
        xfile = "p" // procname // ".grid"
        qfile = "p" // procname // ".T"
        CALL readPlot3D(blocks, xfile, qfile)
    END SUBROUTINE read_config


    SUBROUTINE plot3D(blocks, NBLKS, xfile, qfile)
        ! write plt 2d file given blocks, number of blocks,
        ! x and q file names (no file extension), and the bounds for writing
        ! (0 for real grid, 1 to include ghosts)
        IMPLICIT NONE

        TYPE(BLKTYPE) :: blocks(:)
        INTEGER :: IBLK, I, J, NBLKS, bound = 1
        ! OUTPUT FILES (without file exension)
        CHARACTER(20) :: xfile, qfile

        ! FORMAT STATEMENTS
            ! I --> Integer, number following is number of sig figs
            ! E --> scientific notation,
                        ! before decimal is sig figs of exponent?
                        ! after decimal is sig figs of value
            ! number before letter is how many entries on single line
                ! before newline (number of columns)
        10     FORMAT(I10)
        20     FORMAT(10I10)
        30     FORMAT(10E20.8)

        !!! FORMATTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! OPEN FILES
        OPEN(UNIT=gridUnit,FILE = TRIM(xfile) // '.form.xyz',FORM='formatted')
        OPEN(UNIT=tempUnit,FILE = TRIM(qfile) // '.form.dat',FORM='formatted')

        ! WRITE TO GRID FILE
        WRITE(gridUnit, 10) NBLKS
        WRITE(gridUnit, 20) ( IMAXBLK, JMAXBLK, IBLK=1, NBLKS)
!         WRITE(gridUnit, 20) ( blocks(IBLK)%IMAX, blocks(IBLK)%JMAX, IBLK=1, NBLK)
        DO IBLK = 1, NBLKS
            WRITE(gridUnit, 30) ( (blocks(IBLK)%mesh%x(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound), &
                                ( (blocks(IBLK)%mesh%y(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound)
        END DO


        ! WRITE TO TEMPERATURE FILE
            ! When read in paraview, 'density' will be equivalent to temperature
        WRITE(tempUnit, 10) NBLKS
        WRITE(tempUnit, 20) ( IMAXBLK, JMAXBLK, IBLK=1, NBLKS)
        DO IBLK = 1, NBLKS

            WRITE(tempUnit, 30) tRef,dum,dum,dum
            WRITE(tempUnit, 30) ( (blocks(IBLK)%mesh%T(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound)
        END DO

        ! CLOSE FILES
        CLOSE(gridUnit)
        CLOSE(tempUnit)

        !!! UNFORMATTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! OPEN FILES
        OPEN(UNIT=gridUnit,FILE = TRIM(xfile) // '.xyz',FORM='unformatted')
        OPEN(UNIT=tempUnit,FILE = TRIM(qfile) // '.dat',FORM='unformatted')

        ! WRITE TO GRID FILE (UNFORMATTED)
            ! (Paraview likes unformatted better)
        WRITE(gridUnit) NBLKS
        WRITE(gridUnit) ( IMAXBLK, JMAXBLK, IBLK=1, NBLKS)
!         WRITE(gridUnit) ( blocks(IBLK)%IMAX, blocks(IBLK)%JMAX, IBLK=1, NBLK)
        DO IBLK = 1, NBLKS
            WRITE(gridUnit) ( (blocks(IBLK)%mesh%x(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                            ( (blocks(IBLK)%mesh%y(I,J), I=1,IMAXBLK), J=1,JMAXBLK)
        END DO


        ! WRITE TO TEMPERATURE FILE
            ! When read in paraview, 'density' will be equivalent to temperature
        WRITE(tempUnit) NBLKS
        WRITE(tempUnit) ( IMAXBLK, JMAXBLK, IBLK=1, NBLKS)
        DO IBLK = 1, NBLKS

            WRITE(tempUnit) tRef,dum,dum,dum
            WRITE(tempUnit) ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK)
        END DO

        ! CLOSE FILES
        CLOSE(gridUnit)
        CLOSE(tempUnit)
    END SUBROUTINE plot3D

    SUBROUTINE readPlot3D(blocks, xfile, qfile)
        IMPLICIT NONE

        TYPE(BLKTYPE) :: blocks(:)
        INTEGER :: IBLK, I, J, NBLKS
        INTEGER :: NBLKREAD, IMAXBLKREAD, JMAXBLKREAD, bound = 1
        REAL(KIND=8) :: dum1, dum2, dum3, dum4
        ! OUTPUT FILES (without file exension)
        CHARACTER(20) :: xfile, qfile

        ! FORMAT STATEMENTS
            ! I --> Integer, number following is number of sig figs
            ! E --> scientific notation,
                        ! before decimal is sig figs of exponent?
                        ! after decimal is sig figs of value
            ! number before letter is how many entries on single line
                ! before newline (number of columns)
        10     FORMAT(I10)
        20     FORMAT(10I10)
        30     FORMAT(10E20.8)

        !!! FORMATTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! OPEN FILES
        OPEN(UNIT=gridUnit,FILE = TRIM(xfile) // '.form.xyz',FORM='formatted')
        OPEN(UNIT=tempUnit,FILE = TRIM(qfile) // '.form.dat',FORM='formatted')

        ! READ GRID FILE
        READ(gridUnit, 10) NBLKREAD
        READ(gridUnit, 20) ( IMAXBLKREAD, JMAXBLKREAD, IBLK=1, NBLKREAD)
        DO IBLK = 1, NBLKREAD
            READ(gridUnit, 30) ( (blocks(IBLK)%mesh%x(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound), &
                               ( (blocks(IBLK)%mesh%y(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound)
        END DO

        ! READ TEMPERATURE FILE
        READ(tempUnit, 10) NBLKREAD
        READ(tempUnit, 20) ( IMAXBLKREAD, JMAXBLKREAD, IBLK=1, NBLKREAD)
        DO IBLK = 1, NBLKREAD

!             READ(tempUnit, 30) tRef,dum,dum,dum
            READ(tempUnit, 30) dum1, dum2, dum3, dum4
            READ(tempUnit, 30) ( (blocks(IBLK)%mesh%T(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound), &
                               ( (blocks(IBLK)%mesh%T(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound), &
                               ( (blocks(IBLK)%mesh%T(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound), &
                               ( (blocks(IBLK)%mesh%T(I,J), I=1-bound,IMAXBLK+bound), J=1-bound,JMAXBLK+bound)
        END DO

        ! CLOSE FILES
        CLOSE(gridUnit)
        CLOSE(tempUnit)


    END SUBROUTINE readPlot3D

!     SUBROUTINE readPlot3D(blocks)
!         IMPLICIT NONE

!         TYPE(BLKTYPE) :: blocks(:)
!         INTEGER :: IBLK, I, J
!         ! READ INFO FOR BLOCK DIMENSIONS
!         INTEGER :: NBLKREAD, IMAXBLKREAD, JMAXBLKREAD
!         ! OUTPUT FILES
!         CHARACTER(20) :: xfile, qfile

!         ! FORMAT STATEMENTS
!             ! I --> Integer, number following is number of sig figs
!             ! E --> scientific notation,
!                         ! before decimal is sig figs of exponent?
!                         ! after decimal is sig figs of value
!             ! number before letter is how many entries on single line
!                 ! before newline (number of columns)
!         10     FORMAT(I10)
!         20     FORMAT(10I10)
!         30     FORMAT(10E20.8)

!         !!! FORMATTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         ! OPEN FILES
! !         OPEN(UNIT=gridUnit,FILE= TRIM(casedir) // 'grid_form.xyz',FORM='formatted')
! !         OPEN(UNIT=tempUnit,FILE= TRIM(casedir) // 'T_form.dat',FORM='formatted')
!         OPEN(UNIT=gridUnit,FILE= 'grid_form.xyz',FORM='formatted')
!         OPEN(UNIT=tempUnit,FILE= 'T_form.dat',FORM='formatted')

!         ! READ GRID FILE
!         READ(gridUnit, 10) NBLKREAD
!         READ(gridUnit, 20) ( IMAXBLKREAD, JMAXBLKREAD, IBLK=1, NBLKREAD)
! !         WRITE(gridUnit, 20) ( blocks(IBLK)%IMAX, blocks(IBLK)%JMAX, IBLK=1, NBLK)
!         DO IBLK = 1, NBLKREAD
!             READ(gridUnit, 30) ( (blocks(IBLK)%mesh%x(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
!                                 ( (blocks(IBLK)%mesh%y(I,J), I=1,IMAXBLK), J=1,JMAXBLK)
!         END DO


!         ! READ TEMPERATURE FILE
!             ! When read in paraview, 'density' will be equivalent to temperature
!         READ(tempUnit, 10) NBLKREAD
!         READ(tempUnit, 20) ( IMAXBLKREAD, JMAXBLKREAD, IBLK=1, NBLKREAD)
!         DO IBLK = 1, NBLKREAD

!             READ(tempUnit, 30) tRef,dum,dum,dum
!             READ(tempUnit, 30) ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
!                                 ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
!                                 ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
!                                 ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK)
!         END DO

!         ! CLOSE FILES
!         CLOSE(gridUnit)
!         CLOSE(tempUnit)
!     END SUBROUTINE readPlot3D

    SUBROUTINE compositePlot3D()
        type(blktype), ALLOCATABLE :: blocks(:)
        type(proctype), target :: procs(nprocs)
        type(proctype), pointer :: p
        CHARACTER(2) :: procname
        CHARACTER(20) :: xfile, qfile

        integer :: procsort(NBLK), IDsSort(NBLK), I, ii
        allocate(blocks(NBLK))
        ! read block amalgamation file
        OPEN(UNIT=55,FILE = 'blockrebuild.dat',FORM='formatted')
        read(55,*)
        do I = 1, NBLK
            read(55,*) Ii, procsort(I), IDsSort(I)
        end do
        CLOSE(55)

        OPEN(UNIT=65,FILE = 'procrebuild.dat',FORM='formatted')

        do i = 1, NPROCs
            p => procs(I)
            READ(65,*) p%NBLK
            allocate(p%blocks(p%NBLK))

            IF (p%ID<10) THEN
                ! IF SINGLE DIGIT, PAD WITH 0 IN FRONT
                WRITE(procname, '(A,I1)') '0', p%ID
            ELSE
                WRITE(procname, '(I2)') p%ID
            END IF
            xfile = "p" // procname // ".grid"
            qfile = "p" // procname // ".T"
            call readplot3d(p%blocks, xfile, qfile)

        end do
        CLOSE(65)

        do i = 1, nblk
            blocks(I) = procs( procsort(i) )%blocks( idsSort(i) )
        end do


        call plot3d(blocks, nblk, 'grid', 'T')



    END SUBROUTINE compositePlot3D


    SUBROUTINE write_res(res_hist)
        TYPE(RESLIST), POINTER :: res_hist
        ! pointer to iterate linked list
        TYPE(RESLIST), POINTER :: hist

        ! open residual file
!         OPEN(UNIT=resUnit,FILE= TRIM(casedir) // 'res_hist.dat')
        OPEN(UNIT=resUnit,FILE = 'res_hist.dat')
        ! column headers
        WRITE(resUnit,*) 'ITER      RESID'

        ! point to residual linked list
        hist => res_hist
        ! skip first link, empty from iteration loop design
        hist => hist%next
        ! write residual history to file until list ends
        DO
            IF ( .NOT. ASSOCIATED(hist) ) EXIT
            ! write iteration and residual in two columns
            WRITE(resUnit,*) hist%iter, hist%res
            hist => hist%next
        END DO

        CLOSE(resUnit)
    END SUBROUTINE write_res


END MODULE IO
