       program hello
       implicit none
       integer*4 :: nthreads, tid, omp_get_num_threads, omp_get_thread_num

       !$OMP PARALLEL PRIVATE(NTHREADS,TID)

       tid = omp_get_thread_num()
       print*, 'Hello Dr. D from thread = ', TID

       if(tid.eq.0)then
          nthreads = omp_get_num_threads()
          print*, 'Number of threads = ', nthreads
       end if

       !$OMP END PARALLEL

       end program hello
                 
