program leastSquares

! MAE 267
! HW2
! Problem 1 - Least-Squares Linear Fit
! Logan Halstrom
! 28 September 2015

! DESCRIPTION:  Calculate and plot least-squares linear fit

!   No undeclared variables
    implicit none

!   DECLARE VARIABLES
!   Length of Input Vectors
    real, parameter :: n = 20
!   Pendulum Parameters
    real, dimensions(1:n) :: x, y

!   GIVEN INPUT DATA
    x = (/ -4.91, -3.84, -2.41, -2.62, -3.78, -0.52, -1.83, &
    -2.01, 0.28, 1.08, -0.94, 0.59, 0.69, 3.04, 1.01, 3.60, &
    4.53, 5.13, 4.43, 4.12 /)

    y = (/ -8.18, -7.49, -7.11, -6.15, -5.62, -3.30, -2.05, &
    -2.83, -1.16, 0.52, 0.21, 1.73, 3.96, 4.26, 5.75,     &
    6.67, 7.70, 7.31, 9.05, 10.95 /)

end program leastSquares


