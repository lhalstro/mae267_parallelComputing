program hw2

! MAE 267
! HW2
! Problem 1 - Least-Squares Linear Fit
! Logan Halstrom
! 28 September 2015

! DESCRIPTION:  Calculate and plot least-squares linear fit

!   No undeclared variables
    implicit none

!   DECLARE VARIABLES
    integer, parameter :: n = 20, x_unit = 21, y_unit = 22
    integer :: i
!   Pendulum Parameters
    real, dimension(1:n) :: x, y

    interface
        subroutine leastSquares(x, y)
            real, dimension(:) :: x, y
        end subroutine leastSquares
    end interface

!   GIVEN INPUT DATA
    x = (/ -4.91, -3.84, -2.41, -2.62, -3.78, -0.52, -1.83, &
    -2.01, 0.28, 1.08, -0.94, 0.59, 0.69, 3.04, 1.01, 3.60, &
    4.53, 5.13, 4.43, 4.12 /)

    y = (/ -8.18, -7.49, -7.11, -6.15, -5.62, -3.30, -2.05, &
    -2.83, -1.16, 0.52, 0.21, 1.73, 3.96, 4.26, 5.75,     &
    6.67, 7.70, 7.31, 9.05, 10.95 /)

    ! Write input data to file
!     open (unit=x_unit,file="x.dat",action="write",status="replace")
!     write (x_unit,*) x
!     close (x_unit)
!     open (unit=y_unit,file="y.dat",action="write",status="replace")
!     write (y_unit,*) y
!     close (y_unit)

    open (unit=x_unit,file="x.dat",action="write",status="replace")
    do i=1,n
        ! write newline separated file (delimiter is contained in '')
        write (x_unit,"(32(f0.6,'',:))") x(i)
    end do
    close (x_unit)
    open (unit=y_unit,file="y.dat",action="write",status="replace")
    do i=1,n
        write (y_unit,"(32(f0.6,'',:))") y(i)
    end do
    close (y_unit)

    call leastSquares(x,y)

end program hw2

subroutine leastSquares(x, y)
    implicit none

    integer :: nx, ny, out_unit = 20
    real, dimension(:) :: x, y
    real sumx, sumy, sumxx, sumyy, sumxy, avex, avey, m, b, r

    nx = size(x)
    ny = size(y)
    sumx = sum(x)
    sumy = sum(y)
    sumxx = sum(x*x)
    sumyy = sum(y*y)
    sumxy = sum(x*y)
    avex = sum(x)/nx
    avey = sum(y)/ny

    ! Perform Least Squares Fit
    m = (sumxy - sumx * avey) / (sumxx - sumx * avex)
    b = avey - (m * avex)

    ! Correlation Coeff
    r = (nx * sumxy - sumx * sumy) / &
            sqrt((nx * sumxx - sumx ** 2) * (nx * sumyy - sumy ** 2))

    ! Output Results
    ! Create i/o output file
    open (unit=out_unit,file="output.txt",action="write",status="replace")
    write (out_unit,*) 'Least-squares linear fit results:'
    write (out_unit,*) 'm=', m, 'b=', b
    write (out_unit,*) 'Correlation Coefficient:'
    write (out_unit,*) 'r=', r
    close (out_unit)

end subroutine

