      program master
c
c     this program solves the 1d heat conduction equation with a 
c     point-jacobi routine by calling programs code1 and code2 that
c     have been changed to become subroutines. code1 and code2 each
c     use 2 processors.  The global domain extends from 0 to 1 so that
c     each processor has 1 block that covers 1/4 of the overall domain.
c
      implicit none
      include "mpif.h"
c
      real*8, allocatable, dimension(:) :: x,t,dt
      real*8  :: dtmax,dtmax2
      integer :: i,npts,iter,niter,imax,imax2,icode,iproc,iproc2
      integer :: nprocs,myid,mykey
      integer :: nprocs1,myid1,root1,nprocs2,myid2,root2
      integer :: status,ierror
      integer :: MPI_COMM_CODE1,MPI_COMM_CODE2
      integer :: MPI_COMM_CODE12,MPI_COMM_CODE21
      integer :: MPI_GROUP_WORLD
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c     initialize mpi for master
c
      call MPI_Init(IERROR)
      call MPI_Comm_Rank(MPI_COMM_WORLD,MYID,IERROR)
      call MPI_Comm_Size(MPI_COMM_WORLD,NPROCS,IERROR)
      call MPI_Comm_Group(MPI_COMM_WORLD,MPI_GROUP_WORLD,IERROR)
      print*,'mpi initialized by',myid,' for master'
      print*,'myid for world',myid
      print*,'nprocs for world',nprocs
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
c
      if(myid<=1) then
c
c        create the code1 communicator with 2 processors
c
         mykey = 1
         call MPI_Comm_Split(MPI_COMM_WORLD,mykey,myid,
     &                       MPI_COMM_CODE1,ierror)
      else
c
c        create the code2 communicator with remainder of processors
c
         mykey = 2
         call MPI_Comm_Split(MPI_COMM_WORLD,mykey,myid,
     &                       MPI_COMM_CODE2,ierror)
      endif
      if(myid==1.or.myid==2) then
c
c        create the communicator for inter-code communication
c        note that processors 1 and 2 are known to be at the
c        inter-code interface
c
         mykey = 3
         call MPI_Comm_Split(MPI_COMM_WORLD,mykey,myid,
     &                       MPI_COMM_CODE12,ierror)
      endif
c
      if(myid<=1) then
         call MPI_Comm_Rank(MPI_COMM_CODE1,myid1,ierror)
         call MPI_Comm_Size(MPI_COMM_CODE1,nprocs1,ierror)
         if(myid1==0) root1 = myid1
         print*,'myid for code1',myid1
         print*,'nprocs for code1',nprocs1
      else
         call MPI_Comm_Rank(MPI_COMM_CODE2,myid2,ierror)
         call MPI_Comm_Size(MPI_COMM_CODE2,nprocs2,ierror)
         if(myid2==0) root2 = myid2
         print*,'myid for code2',myid2
         print*,'nprocs for code2',nprocs2
      endif
c
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      stop
c
      call MPI_Finalize(IERROR)
c
      end program
