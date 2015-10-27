! MAE 267
! LOGAN HALSTROM
! 12 OCTOBER 2015

! DESCRIPTION:  This module creates a grid and temperature file in
!               the plot3D format for steady state solution

MODULE plot3D_module
    USE CONSTANTS
    USE BLOCKMOD
    IMPLICIT NONE

    ! VARIABLES
    INTEGER :: gridUnit  = 30   ! Unit for grid file
    INTEGER :: tempUnit = 21    ! Unit for temp file
    REAL(KIND=8) :: tRef = 1.D0          ! tRef number
    REAL(KIND=8) :: dum = 0.D0          ! dummy values
    INTEGER :: nBlocks = 1      ! number of blocks

    CONTAINS
    SUBROUTINE plot3D(mesh)
        IMPLICIT NONE

        TYPE(MESHTYPE) :: mesh
        INTEGER :: i, j

        ! FORMAT STATEMENTS
        10     FORMAT(I10)
        20     FORMAT(10I10)
        30     FORMAT(10E20.8)

!         ! OPEN FILES
!         OPEN(UNIT=gridUnit,FILE='grid.xyz',FORM='formatted')
!         OPEN(UNIT=tempUnit,FILE='temperature.dat',FORM='formatted')

!         ! WRITE TO GRID FILE (FORMATTED)
!         WRITE(gridUnit,10) nBlocks
!         WRITE(gridUnit,20) IMAX,JMAX
!         WRITE(gridUnit,30) ((mesh(i,j)%x,i=1,IMAX),j=1,JMAX), ((mesh(i,j)%y,i=1,IMAX),j=1,JMAX)

!         ! WRITE TO TEMPERATURE FILE
!         WRITE(tempUnit,10) nBlocks
!         WRITE(tempUnit,20) IMAX,JMAX
!         WRITE(tempUnit,30) tRef,dum,dum,dum
!         WRITE(tempUnit,30) ((mesh(i,j)%T,i=1,IMAX),j=1,JMAX), ((mesh(i,j)%T,i=1,IMAX),j=1,JMAX), &
!                            ((mesh(i,j)%T,i=1,IMAX),j=1,JMAX), ((mesh(i,j)%T,i=1,IMAX),j=1,JMAX)

        ! OPEN FILES
        OPEN(UNIT=gridUnit,FILE='grid.xyz',FORM='unformatted')
        OPEN(UNIT=tempUnit,FILE='temperature.dat',FORM='unformatted')

        ! WRITE TO GRID FILE (UNFORMATTED)
            ! (Paraview likes unformatted better)
        WRITE(gridUnit) nBlocks
        WRITE(gridUnit) IMAX,JMAX
        WRITE(gridUnit) ((mesh%x(i,j),i=1,IMAX),j=1,JMAX), ((mesh%y(i,j),i=1,IMAX),j=1,JMAX)

        ! WRITE TO TEMPERATURE FILE
            ! When read in paraview, 'density' will be equivalent to temperature
        WRITE(tempUnit) nBlocks
        WRITE(tempUnit) IMAX,JMAX
        WRITE(tempUnit) tRef,dum,dum,dum
        WRITE(tempUnit) ((mesh%T(i,j),i=1,IMAX),j=1,JMAX), ((mesh%T(i,j),i=1,IMAX),j=1,JMAX), &
                           ((mesh%T(i,j),i=1,IMAX),j=1,JMAX), ((mesh%T(i,j),i=1,IMAX),j=1,JMAX)

        ! CLOSE FILES
        CLOSE(gridUnit)
        CLOSE(tempUnit)
    END SUBROUTINE plot3D
END MODULE plot3D_module
