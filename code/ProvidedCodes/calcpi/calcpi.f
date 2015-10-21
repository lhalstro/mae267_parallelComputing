      PROGRAM CALCPI
      IMPLICIT NONE
      REAL(kind=8) :: A,B,PI
      INTEGER :: N
      A = 0.0d0   !DEFINE INTERVAL START
      B = 1.0d0   !DEFINE INTERVAL STOP
      !READ THE NUMBER OF SUB-INTERVALS
      PRINT *,'INPUT THE NUMBER OF SUB-INTERVALS'
      READ(*,*) N
        CALL TRAP(A,B,N,PI)
        PRINT *,'PI = ',PI
      STOP
      END
