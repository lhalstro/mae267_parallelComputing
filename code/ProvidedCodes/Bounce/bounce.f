C-------------------------------------------------------------------
C     This program times blocking send/receives, and reports the
C     latency and bandwidth of the communication system.  It is
C     designed to run on an even number of nodes.
C-------------------------------------------------------------------
      program bounce
      parameter (maxbytes=1000000)
      parameter (maxcount=1000000)
      parameter (nRepeats=10)
      parameter (nsizes=7)
      implicit real*8 (a-h,o-z)
      include "mpif.h"
      dimension sbuf(maxcount), rbuf(maxcount)
      dimension length(nsizes)
      integer   status(MPI_STATUS_SIZE)
C---------------------------------------
C     define an array of message lengths
C---------------------------------------
      length(1) = 0
      length(2) = 1
      length(2) = 10
      length(3) = 100
      length(4) = 1000
      length(5) = 10000
      length(6) = 100000
      length(7) = 1000000
C---------------------------------------
C     initialize the send buffer to zero
C---------------------------------------
      do n=1,maxcount
      sbuf(n) = 0.0D0
      rbuf(n) = 0.0D0
      end do
C-------------------------------------
C     set up the parallel environment
C-------------------------------------
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world,nNodes,ierr)
      call mpi_comm_rank(mpi_comm_world,nodeID,ierr)
      if (nodeID.eq.0) write(*,*)'Number of processors =',nNodes
      call system('hostname')
C
C     print relevant environment variables
C
      if (mod(nNodes,2) .ne. 0) then
         if (nodeID .eq. 0) then
            write(6,*) ' You must specify an even number of nodes.'
         end if
         call mpi_finalize(ierr)
      end if
C-------------------------------------------------------------
C     send/recv to ensure that the routines are loaded
C-------------------------------------------------------------
      ns = 2
      if (mod(nodeID,2) .eq. 0) then
         call mpi_send(sbuf, length(ns), MPI_REAL8, nodeID+1, 1,
     .                 MPI_COMM_WORLD, ierr)
      else
         call mpi_recv(rbuf, length(ns), MPI_REAL8, nodeID-1, 1, 
     .                 MPI_COMM_WORLD, status, ierr)
      end if
      if (mod(nodeID,2) .eq. 1) then
         call mpi_send(sbuf, length(ns), MPI_REAL8, nodeID-1, 1,
     .                 MPI_COMM_WORLD, ierr)
      else
         call mpi_recv(rbuf, length(ns), MPI_REAL8, nodeID+1, 1,
     .                 MPI_COMM_WORLD, status, ierr)
      end if
C---------------------------------------------------------
C     send or receive messages, and time it.
C     even nodes send, odd nodes receive, then the reverse
C---------------------------------------------------------
      do ns=1, nsizes
         time1 = MPI_WTIME()
         do nr=1, nRepeats
C----------------------------------------------
C           send in one direction i->i+1
C----------------------------------------------
         if (mod(nodeID,2) .eq. 0) then
         call mpi_send(sbuf, length(ns), MPI_REAL8, nodeID+1, 1,
     .                 MPI_COMM_WORLD, ierr)
         else
         call mpi_recv(rbuf, length(ns), MPI_REAL8, nodeID-1, 1,
     .                 MPI_COMM_WORLD, status, ierr)
         end if
C---------------------------------------------------
C           send in the reverse direction i+1->i
C---------------------------------------------------
         if (mod(nodeID,2) .eq. 1) then
         call mpi_send(sbuf, length(ns), MPI_REAL8, nodeID-1, 1,
     .                 MPI_COMM_WORLD, ierr)
         else
         call mpi_recv(rbuf, length(ns), MPI_REAL8, nodeID+1, 1,
     .                 MPI_COMM_WORLD, status, ierr)
         end if
         end do
         time2 = MPI_WTIME()
         if (nodeID .eq. 0) then
         write(6,fmt='(A,I9,A,10X,A,F10.4,A)') 'msglen =',8*length(ns),
     &   ' bytes,','elapsed time =',0.5D3*(time2-time1)/nRepeats,' msec'
         call flush(6)
         end if
         if (ns .eq. 1) then
            tlatency = 0.5D6*(time2-time1)/nRepeats
         end if
         if (ns .eq. nsizes) then
            bw = 8.*length(ns)/(0.5D6*(time2-time1)/nRepeats)
         end if
      end do
C---------------------------------------------------------
C     report apporximate numbers for bandwidth and latency
C---------------------------------------------------------
      if (nodeID .eq. 0) then
         write(6,fmt='(A,F6.1,A)') 'latency =',tlatency,' microseconds'
         write(6,*) 'bandwidth =',bw,' MBytes/sec'
c        write(6,fmt='(A,F8.4,A)') 'bandwidth =',bw,' MBytes/sec'
         write(6,fmt='(A)') '(approximate values for mp_bsend/mp_brecv)'
      end if
 11   continue
      call mpi_finalize(ierr)
      end
