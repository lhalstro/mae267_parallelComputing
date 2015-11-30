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

    SUBROUTINE read_config(b)
        ! Read block connectivity file with neighbor and BC info
        ! for each processor.
        ! Also read PLOT3D restart files for each processor.

!         TYPE(PROCTYPE), TARGET :: procs(:)
!         TYPE(PROCTYPE), POINTER :: p
!         ! BLOCK DATA TYPE
!         TYPE(BLKTYPE), POINTER :: b
!         INTEGER :: IP, IB, BLKFILE = 99
!         CHARACTER(2) :: procname
!         CHARACTER(20) :: xfile, qfile

!         33 FORMAT(A)
!         11 FORMAT( 3I7)
!         22 FORMAT(33I7)
!         44 FORMAT(33A7)





        ! BLOCK DATA TYPE
        TYPE(BLKTYPE) :: b(:)
        INTEGER :: I, BLKFILE = 99
        ! READ INFOR FOR BLOCK DIMENSIONS
        INTEGER :: NBLKREAD, IMAXBLKREAD, JMAXBLKREAD

        11 FORMAT(3I5)
        33 FORMAT(A)
        22 FORMAT(33I5)
        44 FORMAT(33A5)

        OPEN (UNIT = BLKFILE , FILE = "blockconfig.dat", form='formatted')
        ! WRITE AMOUNT OF BLOCKS AND DIMENSIONS
        READ(BLKFILE,*)
        READ(BLKFILE, 11) NBLK, IMAXBLK, JMAXBLK
        READ(BLKFILE,*)
        DO I = 1, NBLK
            ! FOR EACH BLOCK, WRITE BLOCK NUMBER, STARTING/ENDING GLOBAL INDICES.
            ! THEN BOUNDARY CONDITION AND NEIGHBOR NUMBER FOR EACH FACE:
            ! NORTH EAST SOUTH WEST
            READ(BLKFILE, 22) b(I)%ID, &
                b(I)%IMIN, b(I)%JMIN, &
                b(I)%NB%N, &
                b(I)%NB%NE, &
                b(I)%NB%E, &
                b(I)%NB%SE, &
                b(I)%NB%S, &
                b(I)%NB%SW, &
                b(I)%NB%W, &
                b(I)%NB%NW, &
                b(I)%ORIENT
        END DO
        CLOSE(BLKFILE)
    END SUBROUTINE read_config


    SUBROUTINE plot3D(blocks, NBLKS, xfile, qfile)
        IMPLICIT NONE

        TYPE(BLKTYPE) :: blocks(:)
        INTEGER :: IBLK, I, J, NBLKS
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
            WRITE(gridUnit, 30) ( (blocks(IBLK)%mesh%x(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%y(I,J), I=1,IMAXBLK), J=1,JMAXBLK)
        END DO


        ! WRITE TO TEMPERATURE FILE
            ! When read in paraview, 'density' will be equivalent to temperature
        WRITE(tempUnit, 10) NBLKS
        WRITE(tempUnit, 20) ( IMAXBLK, JMAXBLK, IBLK=1, NBLKS)
        DO IBLK = 1, NBLKS

            WRITE(tempUnit, 30) tRef,dum,dum,dum
            WRITE(tempUnit, 30) ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK)
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

    SUBROUTINE readPlot3D(blocks)
        IMPLICIT NONE

        TYPE(BLKTYPE) :: blocks(:)
        INTEGER :: IBLK, I, J
        ! READ INFO FOR BLOCK DIMENSIONS
        INTEGER :: NBLKREAD, IMAXBLKREAD, JMAXBLKREAD
        ! OUTPUT FILES
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
!         OPEN(UNIT=gridUnit,FILE= TRIM(casedir) // 'grid_form.xyz',FORM='formatted')
!         OPEN(UNIT=tempUnit,FILE= TRIM(casedir) // 'T_form.dat',FORM='formatted')
        OPEN(UNIT=gridUnit,FILE= 'grid_form.xyz',FORM='formatted')
        OPEN(UNIT=tempUnit,FILE= 'T_form.dat',FORM='formatted')

        ! READ GRID FILE
        READ(gridUnit, 10) NBLKREAD
        READ(gridUnit, 20) ( IMAXBLKREAD, JMAXBLKREAD, IBLK=1, NBLKREAD)
!         WRITE(gridUnit, 20) ( blocks(IBLK)%IMAX, blocks(IBLK)%JMAX, IBLK=1, NBLK)
        DO IBLK = 1, NBLKREAD
            READ(gridUnit, 30) ( (blocks(IBLK)%mesh%x(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%y(I,J), I=1,IMAXBLK), J=1,JMAXBLK)
        END DO


        ! READ TEMPERATURE FILE
            ! When read in paraview, 'density' will be equivalent to temperature
        READ(tempUnit, 10) NBLKREAD
        READ(tempUnit, 20) ( IMAXBLKREAD, JMAXBLKREAD, IBLK=1, NBLKREAD)
        DO IBLK = 1, NBLKREAD

            READ(tempUnit, 30) tRef,dum,dum,dum
            READ(tempUnit, 30) ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK), &
                                ( (blocks(IBLK)%mesh%T(I,J), I=1,IMAXBLK), J=1,JMAXBLK)
        END DO

        ! CLOSE FILES
        CLOSE(gridUnit)
        CLOSE(tempUnit)
    END SUBROUTINE readPlot3D



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
