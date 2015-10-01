program pendulumPeriod

! MAE 267
! HW1
! Problem 1 -
! Logan Halstrom
! 28 September 2015

! DESCRIPTION:  Calculate period of a pendulum

!   No undeclared variables
    implicit none

!   DECLARE VARIABLES
!   Pendulum Parameters
    real :: L, T
!   Constants
    real, parameter :: pi = 3.141592653589793, g = 9.81
!   Input/Output
    character(1) :: yesno
!   Number for i/o output
    integer, parameter :: out_unit=20

!   Create i/o output file
    open (unit=out_unit,file="output.txt",action="write",status="replace")

    main_loop: do
!       GET PENDULUM LENGTH FROM USER
        write (*,*) 'Enter pendulum length (m):'
        write (out_unit,*) 'Enter pendulum length (m):'
        read (*,*) L
        write (out_unit,*) L

!       CALCULATE PENDULUM PERIOD
        T = 2 * pi * sqrt(L / g)

        write (*,*) 'Period of a pendulum of L=', L, '(m) is T=', T, 's'
        write (out_unit,*) 'Period of a pendulum of L=', L, '(m) is T=', T, 's'

!       PERFORM CALCULATIONS AGAIN IF DESIRED
        write (*,*) 'Perform more calculations? (y/n)'
        write (out_unit,*) 'Perform more calculations? (y/n)'
        read (*,*) yesno
        write (out_unit,*) yesno
        if (yesno=='n' .or. yesno == 'N') exit main_loop

    end do main_loop
    close (out_unit)
end program pendulumPeriod


