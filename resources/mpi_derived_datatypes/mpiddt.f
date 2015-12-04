      PROGRAM MPIDDT
!
!     THIS PROGRAM IS A TEST OF MPI DERIVED DATATYPE OPERATIONS
!
      IMPLICIT NONE

      include "mpif.h"

      REAL(KIND=8) :: A(10,10),B(10,10),C(5,6)
      INTEGER :: I,J,COLUMNTYPE,SUBBLOCK,IERROR,MYID,NPROCS,SIZEOFREAL
      INTEGER :: STATUS(MPI_STATUS_SIZE)

!     INITIALIZE MPI
      CALL MPI_Init(IERROR)
      CALL MPI_Comm_Rank(MPI_COMM_WORLD,MYID,IERROR)
      PRINT *,'MYID = ',MYID
      CALL MPI_Comm_Size(MPI_COMM_WORLD,NPROCS,IERROR)

      IF(NPROCS/=2) STOP

      A=0.0D0
      B=0.0D0
      C=-1.0D0

!     INITIALIZE A AND B
      IF(MYID==0) THEN
        DO J = 1,10
          DO I = 1,10
            A(I,J) = REAL(I*J)
          ENDDO
        ENDDO
      ELSE
        DO J = 1,10
          DO I = 1,10
            B(I,J) = REAL(I)
          ENDDO
        ENDDO
      ENDIF

!     LOOK AT APPENDIX E P350 OF "USING MPI" BY GROPP ET AL

      CALL MPI_Type_Vector(6,5,10,MPI_DOUBLE_PRECISION,
     &                     SUBBLOCK,IERROR)
      CALL MPI_Type_Commit(SUBBLOCK,IERROR)

      IF(MYID/=0) THEN
        CALL MPI_Send(B(1,1),1,SUBBLOCK,0,0,MPI_COMM_WORLD,
     &                IERROR)
      ELSE
        CALL MPI_Recv(C(1,1),5*6,MPI_DOUBLE_PRECISION,1,0,
     &                MPI_COMM_WORLD,STATUS,IERROR)
      ENDIF

      IF(MYID==0) THEN
        PRINT *,'WRITING OUT C',MYID
        DO J = 6,1,-1
          WRITE(*,10) (C(I,J),I=1,5)
 10       FORMAT(5F10.5)
        ENDDO
      ENDIF

      CALL MPI_FINALIZE(IERROR)

      STOP
      END






