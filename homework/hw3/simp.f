      SUBROUTINE SIMP(A,B,N,INTEGRAL)
      !CALCULATE INTEGRAL OF FUNCTION WITH EVEN INTERVAL (SIMPSONS RULE)
      IMPLICIT NONE
      REAL(kind=8)  :: A,B,F,H,INTEGRAL,X
      INTEGER :: N,I
      ! FUNCTION DEFINITION
      F(X) = 4.0D0/(1.D0+X**2)
      ! INTERVAL DEFINITION (dx)
         ! constant interval
      H = (B-A)/REAL(N)
      ! INITIALIZE INTEGRAL
         ! Loop goes through pairs of odd and even subscript terms
         ! Here, sum the ends of the interval, with the first odd term (x1)
         ! which wouldn't fit into the structure of the loop (a + h) is 2nd term
      INTEGRAL = F(A) + 4.0D0 * F(A + H) + F(B)
      ! LOOP THROUGH SIMPOSON RULE PAIRS
         ! loop though even numbers (2,4,6) and calc odd within loop
      DO I = 2, N-2, 2
         ! start at 3rd term (x2=A+2*dx)
         !CURRENT EVEN TERM
         X = A + REAL(I)*H
         ! add even/odd pairs for each loop iteration (i.e. 2*x2 + 4*x3, etc)
         !x is current even term, x+h is adjacent odd term
         INTEGRAL = INTEGRAL + 2.0D0 * F(X) + 4.0D0 * F(X + H)
      END DO
      INTEGRAL = H / 3.0D0 * INTEGRAL
      RETURN
      END SUBROUTINE SIMP
