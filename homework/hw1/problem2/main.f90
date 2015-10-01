program cableTension

! MAE 267
! HW1
! Problem 2 - Tension of a Cable
! Logan Halstrom
! 28 September 2015

! DESCRIPTION:

!   No undeclared variables
    implicit none

!   DECLARE VARIABLES
!   Pendulum Parameters
    real :: lp=8., lc=8., W=200., dmin=1., dmax=7., dd=0.1
!   Stuff that needs more prescision
    double precision :: T, Tmin=999., d, dopt
!   Number for i/o output
    integer, parameter :: out_unit=20

!   Create i/o output file
    open (unit=out_unit,file="output.txt",action="write",status="replace")

    d = dmin
    do while(d < dmax)
!       TENSION IN CABLE
        T = W * lc * lp / (d * sqrt(lp ** 2 - d ** 2))
        if (T < Tmin) then
!           Minimum Tension
            Tmin = T
!           Length Giving Minimum Tension
            dopt = d
        end if
!       Increment d
        d = d + dd
    end do

    write (*,*) 'Minimum Tension:', Tmin, ' lbs.'
    write (out_unit,*) 'Minimum Tension:', Tmin, ' lbs.'
    write (*,*) 'at optimal distance d=', dopt, ' ft'
    write (out_unit,*) 'at optimal distance d=', dopt, ' ft'

    close (out_unit)

end program cableTension


