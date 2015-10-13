! MAE 267

MODULE CONSTANTS
    IMPLICIT NONE
    ! Pi, grid rotation angle (30 deg)
    REAL(KIND=8), PARAMETER :: pi = 3.141592654D0

CONTAINS
    SUBROUTINE PrintPi(n)
        ! Set size of grid (square)
        REAL(KIND=8) :: n, product
        product = n * pi
        WRITE(*,*) n, 'times pi =', product
    END SUBROUTINE PrintPi
END MODULE CONSTANTS


