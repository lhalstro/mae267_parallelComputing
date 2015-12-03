! MAE 267
! PROJECT 5
! LOGAN HALSTROM
! 29 NOVEMBER 2015


! DESCRIPTION:  Solve heat conduction equation for single block of steel.

! INPUTS: Set grid size, block decomposition, debug in 'config.in'
!         Set number of processors in 'run.sh'

! TO COMPILE:
    ! mpif90 -o main -O3 modules.f90 inout.f90 subroutines.f90 main.f90
        ! makes executable file 'main'
        ! 'rm *.mod' afterward to clean up unneeded compiled files
! TO RUN:
    ! on hpc1 nodes: sbatch run.sh
    ! on hpc1 front end: ./main or ./run.sh



PROGRAM heatTrans
!     USE CLOCK
    USE CONSTANTS
    USE subroutines
    USE IO

    IMPLICIT NONE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! INITIALIZE VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ALL BLOCKS IN ONE LIST
    TYPE(BLKTYPE), ALLOCATABLE :: allblocks(:)

    ! PROCESSORS
    TYPE(PROCTYPE), ALLOCATABLE :: procs(:)
    CHARACTER(2) :: procname
    CHARACTER(20) :: xfile, qfile
    ! ITERATION PARAMETERS
    ! Residual history linked list
    TYPE(RESLIST), POINTER :: res_hist
    ! Maximum number of iterations
    INTEGER :: iter = 1, IBLK, IP
    REAL(KIND=8) :: start_total, end_total
    REAL(KIND=8) :: start_solve, end_solve
    ! CLOCK TOTAL TIME OF RUN
    start_total = MPI_Wtime()


    write(*,*) 'starting mpi'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! START MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! INITIALIZE MPI
    CALL MPI_Init(IERROR)
    ! DETERMINE MY PROCESSOR ID
    ! ARGUMENTS: COMM, MYID, IERROR
    CALL MPI_Comm_rank(MPI_COMM_WORLD, MYID, IERROR)
    ! FIND OUT HOW MANY PROCESSORS ARE USED
    ! ARGUMENTS: COMM, NPROCS, IERROR
    CALL MPI_Comm_size(MPI_COMM_WORLD, NPROCS, IERROR)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! INITIALIZE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! READ INPUTS FROM FILE
    CALL read_input()

    ! have the first processor only set up problem
    IF(MYID == 0) THEN

        write(*,*) 'initializing'

!         ! READ INPUTS FROM FILE
!         CALL read_input()
        ALLOCATE( allblocks(NBLK) )
        ALLOCATE( procs(NPROCS) )
        ! INIITIALIZE GRID SYSTEM
        WRITE(*,*) 'Making mesh...'
        CALL init_gridsystem(allblocks, procs)

        ! CLEAN UP INITIALIZATION
        DEALLOCATE(allblocks, procs)
    END IF

!     ! ONLY PROC 0 READS IN CONFIG DATA, SO BRODCAST TO ALL PROCS
!     ! (syntax: variable to brodcast, size, type, which proc to bcast from, otherstuff)
!     CALL MPI_Bcast(IMAX, 1, MPI_INT, 0, mpi_comm_world, ierror)
!     CALL MPI_Bcast(JMAX, 1, MPI_INT, 0, mpi_comm_world, ierror)

    ! HOLD ALL PROCESSORS UNTIL INITIALIZATION IS COMPLETE
    CALL MPI_Barrier(MPI_COMM_WORLD, IERROR)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SOLVER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! INITIALIZE SOLUTION
    write(*,*) "Initialize for proc ", MYID
    CALL init_solution(blocks, nbrlists, mpilists)

!     if (nprocs == 4) then
!         if (myid == 3) then
!             write(*,*) "block ",   blocks(4)%ID
! !             write(*,*) "iminloc ", blocks(3)%IMINLOC
! !             write(*,*) "Imaxloc ", blocks(3)%IMaxLOC
! !             write(*,*) "jminloc ", blocks(3)%jMINLOC
! !             write(*,*) "jmaxloc ", blocks(3)%jmaxLOC

!             write(*,*) blocks(4)%mesh%term( imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(4)%mesh%V2nd( imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(4)%mesh%V   ( imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(4)%mesh%dt(   imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(4)%mesh%xp(   imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(4)%mesh%x(    imaxblk+1, jmaxblk+1)


!         end if


!     else if (nprocs == 1) then
!         if (myid == 0) then
!             write(*,*) "block ",   blocks(9)%ID
! !             write(*,*) "iminloc ", blocks(10)%IMINLOC
! !             write(*,*) "Imaxloc ", blocks(10)%IMaxLOC
! !             write(*,*) "jminloc ", blocks(10)%jMINLOC
! !             write(*,*) "jmaxloc ", blocks(10)%jmaxLOC
!             write(*,*) blocks(14)%mesh%term( imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(14)%mesh%V2nd( imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(14)%mesh%V   ( imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(14)%mesh%dt(   imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(14)%mesh%xp(   imaxblk+1, jmaxblk+1)
!             write(*,*) blocks(14)%mesh%x(    imaxblk+1, jmaxblk+1)
!             write(*,*) Imax
!         end if

!     end if


    CALL MPI_Barrier(MPI_COMM_WORLD, IERROR)
    ! SOLVE
    WRITE(*,*) 'Solving heat conduction with Processor ', MYID
    CALL solve(blocks, nbrlists, mpilists, iter, res_hist)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SAVE RESULTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    WRITE(*,*) 'Writing results...'

    !TURN THIS ON FOR PJ5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! SAVE SOLUTION AS PLOT3D FILES
        ! MAKE FILE NAME
        IF (MYID<10) THEN
            ! IF SINGLE DIGIT, PAD WITH 0 IN FRONT
            WRITE(procname, '(A,I1)') '0', MYID
        ELSE
            WRITE(procname, '(I2)') MYID
        END IF
        xfile = "p" // procname // ".grid"
        qfile = "p" // procname // ".T"
        CALL plot3D(blocks, MYNBLK, xfile, qfile)
    ! CALC TOTAL WALL TIME
!     end_total = MPI_Wtime()
!     wall_time_total = end_total - start_total

    !TURN THIS ON FOR PJ5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (MYID == 0) THEN
        ! SAVE RESIDUAL HISTORY
        CALL write_res(res_hist)
    END IF
    ! SAVE SOLVER PERFORMANCE PARAMETERS
    CALL output(blocks, iter)


!     if (myid == 0) then
!         call compositePlot3D()
!     end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! CLEAN UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DEALLOCATE(blocks)

    IF (MYID == 0) THEN
        WRITE(*,*) 'Done!'
    END IF

    CALL MPI_Finalize(ierror)

END PROGRAM heatTrans
