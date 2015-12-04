      subroutine code1
c
c     this program solves the 1d heat conduction equation with a 
c     point-jacobi routine and interacts with a similar code called, 
c     code2, to simulate an mpmd environment.  code1 and code2 each
c     use 2 processors.  The global domain extends from 0 to 1 so that
c     each processor has 1 block that covers 1/4 of the overall domain.
c
c     NOTE THAT THIS CODE HAS NOT YET BEEN DEBUGED AND MADE OPERATIONAL
c       AND IS MEANT AS AN EXAMPLE
c
      implicit none
      include "mpif.h"
c
      real*8, allocatable, dimension(:) :: x,t,dt
      real*8  :: dtmax,dtmax2
      integer :: i,npts,iter,niter,imax,imax2,icode,iproc,iproc2
      integer :: nprocs,myid,mykey,mygroup,status,ierror
      integer :: MPI_COMM_CODE1,MPI_COMM_CODE2,MPI_COMM_CODE12
      integer :: MPI_GROUP_WORLD
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c     initialize mpi for code1
c
      call MPI_Init(IERROR)
      call MPI_Comm_Rank(MPI_COMM_WORLD,MYID,IERROR)
      call MPI_Comm_Size(MPI_COMM_WORLD,NPROCS,IERROR)
      call MPI_Comm_Group(MPI_COMM_WORLD,MPI_GROUP_WORLD,IERROR)
      print*,'mpi initialized by',myid,' from code1'
      stop
c
c     create the code1 communicator
c
      mykey = 1
      call MPI_Comm_Split(MPI_COMM_WORLD,mykey,myid,
     &                    MPI_COMM_CODE1,ierror)
      call MPI_Comm_Rank(MPI_COMM_CODE1,myid,ierror)
      call MPI_Comm_Size(MPI_COMM_CODE1,nprocs,ierror)
      call MPI_Comm_Group(MPI_COMM_CODE1,mygroup,ierror)
c
c     send code1 communicator to code2
c
      call MPI_Send(MPI_COMM_CODE1,1,MPI_INTEGER,12,
     &              MPI_COMM_WORLD,ierror)
c
c     receive code2 communicator from code2
c
      call MPI_Recv(MPI_COMM_CODE2,1,MPI_INTEGER,21,
     &              MPI_COMM_WORLD,status,ierror)
c
c     create inter-communicator with code2
c
      call MPI_Intercomm_Create(MPI_COMM_CODE1,0,MPI_COMM_CODE2,0,
     &                          12,MPI_COMM_CODE12,IERROR)
c
      if(myid==0) then
         print*,'inter-communicator successfully created in code1'
      endif
c
      npts = 51
      niter= 2
c
c     initialize grid and temperature (note that ghost nodes are
c       used to deal with block boundary communication)
c
      allocate(x(0:npts+1),t(0:npts+1),dt(0:npts+1))
c
      if(myid==0) then
c        proc 0 of code1 takes care of 0 <= x >= 0.25
         do i=0,npts+1
            x(i) = 0.25*real(i-1)/real(npts-1)
            t(i) = 10.
         enddo
         t(0) = 1.0
      elseif(myid==1) then
c        proc 1 of code1 takes care of 0.25 <= x >= 0.5
         do i=0,npts+1
            x(i) = 0.25+0.25*real(i-1)/real(npts-1)
            t(i) = 10.
         enddo
      endif
c
      do iter = 1,niter
         if(myid==0) then
c
c           message pass information at code1 inter-processor boundary
c           processor 0 of code1 sends to processor 1 of code1
            call MPI_Send(t(npts-1),1,MPI_DOUBLE_PRECISION,0,1,
     &                    MPI_COMM_CODE1,ierror)
c           processor 0 of code1 receives from processor 1 of code1
            call MPI_Recv(t(npts+1),1,MPI_DOUBLE_PRECISION,0,2,
     &                    MPI_COMM_CODE1,status,ierror)
         elseif(myid==1) then
c
c           message pass information at code1 inter-processor boundary
c           processor 1 of code1 sends to processor 0 of code1
            call MPI_Send(t(2),1,MPI_DOUBLE_PRECISION,0,2,
     &                    MPI_COMM_CODE1,ierror)
c           processor 1 of code1 receives from processor 0 of code1
            call MPI_Recv(t(0),1,MPI_DOUBLE_PRECISION,0,1,
     &                    MPI_COMM_CODE1,status,ierror)
c
c           message pass information at code12 inter-processor boundary
c           processor 1 of code1 sends to processor 0 of code2
            call MPI_Send(t(npts-1),1,MPI_DOUBLE_PRECISION,0,5,
     &                    MPI_COMM_CODE12,ierror)
c           processor 1 of code1 receives from processor 0 of code2
            call MPI_Recv(t(npts+1),1,MPI_DOUBLE_PRECISION,0,6,
     &                    MPI_COMM_CODE12,status,ierror)
         endif
c
         do i = 1,npts
            dt(i) = 0.5*(t(i+1)+t(i))
         enddo
c
c        physical boundary conditions
c
         if(myid==0) then
            dt(1) = 0.0d0
         endif
c
c        update temperature
c
         do i = 1,npts
            t(i) = t(i)+dt(i)
         enddo
c
c        check on convergence
c
         dtmax = -1.d10
         do i = 1,npts
            if(abs(dt(i))>dtmax) then
               dtmax = abs(dt(i))
               imax  = i
               iproc = myid
               icode = 1
            endif
         enddo
c
c        receive the dtmax and imax from processor 1 of code1
c
         if(myid==0) then
            call MPI_Recv(dtmax2,1,MPI_DOUBLE_PRECISION,0,90,
     &                       MPI_COMM_CODE1,status,ierror)
            call MPI_Recv(imax2,1,MPI_INTEGER,0,91,
     &                       MPI_COMM_CODE1,status,ierror)
            call MPI_Recv(iproc2,1,MPI_INTEGER,0,92,
     &                       MPI_COMM_CODE1,status,ierror)
            if(dtmax2>dtmax) then
               dtmax = dtmax2
               imax  = imax2
               iproc = iproc2
            endif
         elseif(myid==1) then
            call MPI_Send(dtmax,1,MPI_DOUBLE_PRECISION,0,90,
     &                    MPI_COMM_CODE1,IERROR)
            call MPI_Send(imax,1,MPI_INTEGER,0,91,
     &                    MPI_COMM_CODE1,IERROR)
            call MPI_Send(iproc,1,MPI_INTEGER,0,92,
     &                    MPI_COMM_CODE1,IERROR)
         endif
c
c        synchronize processors for code 1
         call MPI_Barrier(MPI_COMM_CODE1,IERROR)
c
c        receive the dtmax and imax from code2
c
         if(myid==0) then
            call MPI_Recv(dtmax2,1,MPI_DOUBLE_PRECISION,0,100,
     &                       MPI_COMM_CODE12,status,ierror)
            call MPI_Recv(imax2,1,MPI_INTEGER,0,101,
     &                       MPI_COMM_CODE12,status,ierror)
            call MPI_Recv(iproc2,1,MPI_INTEGER,0,102,
     &                       MPI_COMM_CODE12,status,ierror)
c
            if(dtmax2>dtmax) then
               dtmax = dtmax2
               imax  = imax2
               iproc = iproc2
               icode = 2
            endif
c
c           write out the global convergence information
c
            write(6,10) dtmax,imax,iproc,icode
 10         format('dtmax = ',e15.5,2x,'at ',i2,2x,'in proc ',i2,2x,
     &             'in code ',i2)
         endif
c
c        synchronize codes 1 and 2
c
         call MPI_Barrier(MPI_COMM_CODE12,IERROR)
c
      enddo
c
      end subroutine
      
