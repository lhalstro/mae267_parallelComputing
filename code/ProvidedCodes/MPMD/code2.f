      subroutine code2
c
c     this program solves the 1d heat conduction equation with a 
c     point-jacobi routine and interacts with a similar code called, 
c     code1, to simulate an mpmd environment.  code1 and code2 each
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
c     initialize mpi for code2
c
      call MPI_Init(IERROR)
      call MPI_Comm_Rank(MPI_COMM_WORLD,MYID,IERROR)
      call MPI_Comm_Size(MPI_COMM_WORLD,NPROCS,IERROR)
      call MPI_Comm_Group(MPI_COMM_WORLD,MPI_GROUP_WORLD,IERROR)
      print*,'mpi initialized by',myid,' from code2'
      stop
c
c     create the code2 communicator
c
      mykey = 2
      call MPI_Comm_Split(MPI_COMM_WORLD,mykey,myid,
     &                    MPI_COMM_CODE2,ierror)
      call MPI_Comm_Rank(MPI_COMM_CODE2,myid,ierror)
      call MPI_Comm_Size(MPI_COMM_CODE2,nprocs,ierror)
      call MPI_Comm_Group(MPI_COMM_CODE2,mygroup,ierror)
c
c     send code2 communicator to code1
c
      call MPI_Send(MPI_COMM_CODE2,1,MPI_INTEGER,21,
     &              MPI_COMM_WORLD,ierror)
c
c     receive code1 communicator from code1
c
      call MPI_Recv(MPI_COMM_CODE1,1,MPI_INTEGER,12,
     &              MPI_COMM_WORLD,status,ierror)
c
c     create inter-communicator with code1
c
      call MPI_Intercomm_Create(MPI_COMM_CODE2,0,MPI_COMM_CODE1,0,
     &                          12,MPI_COMM_CODE12,IERROR)
c
      if(myid==0) then
         print*,'inter-communicator successfully created for code2'
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
c        proc 0 of code2 takes care of 0.5 <= x >= 0.75
         do i=0,npts+1
            x(i) = 0.50+0.25*real(i-1)/real(npts-1)
            t(i) = 10.
         enddo
         t(0) = 1.0
      elseif(myid==1) then
c        proc 1 of code2 takes care of 0.75 <= x >= 1.0
         do i=0,npts+1
            x(i) = 0.75+0.25*real(i-1)/real(npts-1)
            t(i) = 10.
         enddo
      endif
c
      do iter = 1,niter
         if(myid==0) then
c
c           message pass information at code2 inter-processor boundary
c           processor 0 of code2 sends to processor 1 of code2
            call MPI_Send(t(npts-1),1,MPI_DOUBLE_PRECISION,0,3,
     &                    MPI_COMM_CODE2,ierror)
c           processor 0 of code2 receives from processor 1 of code2
            call MPI_Recv(t(npts+1),1,MPI_DOUBLE_PRECISION,0,4,
     &                    MPI_COMM_CODE2,status,ierror)
c
c           message pass information at code12 inter-processor boundary
c           processor 0 of code2 sends to processor 1 of code1
            call MPI_Send(t(2),1,MPI_DOUBLE_PRECISION,1,6,
     &                    MPI_COMM_CODE12,ierror)
c           processor 0 of code2 receives from processor 1 of code1
            call MPI_Recv(t(0),1,MPI_DOUBLE_PRECISION,1,5,
     &                    MPI_COMM_CODE12,status,ierror)
         elseif(myid==1) then
c
c           message pass information at code2 inter-processor boundary
c           processor 1 of code2 sends to processor 0 of code2
            call MPI_Send(t(2),1,MPI_DOUBLE_PRECISION,0,4,
     &                    MPI_COMM_CODE2,ierror)
c           processor 1 of code2 receives from processor 0 of code2
            call MPI_Recv(t(0),1,MPI_DOUBLE_PRECISION,0,3,
     &                    MPI_COMM_CODE2,status,ierror)
         endif
c
         do i = 1,npts
            dt(i) = 0.5*(t(i+1)+t(i))
         enddo
c
c        physical boundary conditions
c
         if(myid==1) then
            dt(npts) = 0.0d0
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
c        receive the dtmax and imax from processor 1 of code2
c
         if(myid==0) then
            call MPI_Recv(dtmax2,1,MPI_DOUBLE_PRECISION,0,80,
     &                       MPI_COMM_CODE2,status,ierror)
            call MPI_Recv(imax2,1,MPI_INTEGER,0,81,
     &                       MPI_COMM_CODE2,status,ierror)
            call MPI_Recv(iproc2,1,MPI_INTEGER,0,82,
     &                       MPI_COMM_CODE2,status,ierror)
            if(dtmax2>dtmax) then
               dtmax = dtmax2
               imax  = imax2
               iproc = iproc2
            endif
         elseif(myid==1) then
            call MPI_Send(dtmax,1,MPI_DOUBLE_PRECISION,0,80,
     &                    MPI_COMM_CODE2,IERROR)
            call MPI_Send(imax,1,MPI_INTEGER,0,81,
     &                    MPI_COMM_CODE2,IERROR)
            call MPI_Send(iproc,1,MPI_INTEGER,0,82,
     &                    MPI_COMM_CODE2,IERROR)
         endif
c
c        synchronize processors for code 1
         call MPI_Barrier(MPI_COMM_CODE2,IERROR)
c
c        send the dtmax and imax from code2
c
         if(myid==0) then
            call MPI_Send(dtmax,1,MPI_DOUBLE_PRECISION,0,100,
     &                       MPI_COMM_CODE12,IERROR)
            call MPI_Send(imax,1,MPI_INTEGER,0,101,
     &                       MPI_COMM_CODE12,IERROR)
            call MPI_Send(iproc,1,MPI_INTEGER,0,102,
     &                       MPI_COMM_CODE12,IERROR)
         endif
c
c        synchronize codes 1 and 2
c
         call MPI_Barrier(MPI_COMM_CODE12,IERROR)
c
      enddo
c
      end subroutine
      
