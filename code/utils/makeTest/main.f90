! TEST MAKEFILE AND DEPENDANCIES

PROGRAM makeTest
    USE CONSTANTS
    USE subroutines

    IMPLICIT NONE

    WRITE(*,*) "Call 'hello'"
    CALL hello()
    WRITE(*,*) "Print pi from main:", pi
    CALL hello()
    WRITE(*,*) "Print 3*pi with 'UseMod'"
    CALL UseMod(3)

END PROGRAM makeTest
