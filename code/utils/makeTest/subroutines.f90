! MAE 267

MODULE subroutines
    USE CONSTANTS

    IMPLICIT NONE

CONTAINS
    SUBROUTINE hello()
        WRITE(*,*) 'HELLO WORLD!!!'
    END SUBROUTINE hello

    SUBROUTINE UseMod(n)
        REAL(KIND=8) :: n
        CALL PrintPi(n)
    END SUBROUTINE UseMod

END MODULE subroutines


