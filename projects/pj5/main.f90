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

    ! have the first processor only set up problem
    IF(MYID == 0) THEN

        write(*,*) 'initializing'

        ! READ INPUTS FROM FILE
        CALL read_input()
        ALLOCATE( allblocks(NBLK) )
        ALLOCATE( procs(NPROCS) )
        ! INIITIALIZE GRID SYSTEM
        WRITE(*,*) 'Making mesh...'
        CALL init_gridsystem(allblocks, procs)

        ! CLEAN UP INITIALIZATION
        DEALLOCATE(allblocks, procs)
    END IF

    ! HOLD ALL PROCESSORS UNTIL INITIALIZATION IS COMPLETE
    CALL MPI_Barrier(MPI_COMM_WORLD, IERROR)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SOLVER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! INITIALIZE SOLUTION
    write(*,*) "Initialize for proc ", MYID
    CALL init_solution(blocks, nbrlists, mpilists)
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
!     CALL output(blocks, iter)


!     if (myid == 0) then
!         call compositePlot3D()
!     end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! CLEAN UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     DO IBLK = 1, NBLK
!         DEALLOCATE( blocks(IBLK)%mesh%xp   )
!         DEALLOCATE( blocks(IBLK)%mesh%yp   )
!         DEALLOCATE( blocks(IBLK)%mesh%x    )
!         DEALLOCATE( blocks(IBLK)%mesh%y    )
!         DEALLOCATE( blocks(IBLK)%mesh%T    )
!         DEALLOCATE( blocks(IBLK)%mesh%Ttmp )
!         DEALLOCATE( blocks(IBLK)%mesh%dt   )
!         DEALLOCATE( blocks(IBLK)%mesh%V  )
!         DEALLOCATE( blocks(IBLK)%mesh%V2nd )
!         DEALLOCATE( blocks(IBLK)%mesh%term )
!         DEALLOCATE( blocks(IBLK)%mesh%yPP)
!         DEALLOCATE( blocks(IBLK)%mesh%yNP)
!         DEALLOCATE( blocks(IBLK)%mesh%yNN)
!         DEALLOCATE( blocks(IBLK)%mesh%yPN)
!         DEALLOCATE( blocks(IBLK)%mesh%xNN)
!         DEALLOCATE( blocks(IBLK)%mesh%xPN)
!         DEALLOCATE( blocks(IBLK)%mesh%xPP)
!         DEALLOCATE( blocks(IBLK)%mesh%xNP)
!     END DO

    DEALLOCATE(blocks)

    WRITE(*,*) 'Done!'

    CALL MPI_Finalize(ierror)

END PROGRAM heatTrans
