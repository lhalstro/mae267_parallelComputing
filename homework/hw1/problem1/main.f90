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


    main_loop: do
!       GET PENDULUM LENGTH FROM USER
        write (*,*) 'Enter pendulum length (m):'
        read (*,*) L

!       CALCULATE PENDULUM PERIOD
        T = 2 * pi * sqrt(L / g)

        write (*,*) 'Period of a pendulum of L=', L, 'is T=', T, 's'

!       PERFORM CALCULATIONS AGAIN IF DESIRED
        write (*,*) 'Perform more calculations? (y/n)'
        read (*,*) yesno
        if (yesno=='n' .or. yesno == 'N') exit main_loop

    end do main_loop
end program pendulumPeriod


