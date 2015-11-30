! MAE 267
! PROJECT 4
! LOGAN HALSTROM
! 14 NOVEMBER 2015


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

    ! BLOCKS
    TYPE(BLKTYPE), ALLOCATABLE :: blocks(:)
    ! LINKED LISTS STORING NEIGHBOR INFO
    TYPE(NBRLIST) :: nbrlists
    ! PROCESSORS
    TYPE(PROCTYPE), ALLOCATABLE :: procs(:)
    ! ITERATION PARAMETERS
    ! Residual history linked list
    TYPE(RESLIST), POINTER :: res_hist
    ! Maximum number of iterations
    INTEGER :: iter = 1, IBLK
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

!     ! SET NPROCS MANUALLY NOW FOR DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     NPROCS = 4

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! INITIALIZE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! have the first processor only set up problem
    IF(MYID == 0) THEN

        write(*,*) 'initializing'



        ! READ INPUTS FROM FILE
        CALL read_input()
        ALLOCATE( blocks(NBLK) )
        ALLOCATE( procs(NPROCS) )
        ! INIITIALIZE GRID SYSTEM
        WRITE(*,*) 'Making mesh...'
        CALL init_gridsystem(blocks, procs)

        ! CLEAN UP INITIALIZATION
        DEALLOCATE(blocks, procs)
    END IF

    ! HOLD ALL PROCESSORS UNTIL INITIALIZATION IS COMPLETE
    CALL MPI_Barrier(MPI_COMM_WORLD, IERROR)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SOLVER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ! INITIALIZE SOLUTION
!     CALL init_solution(blocks, nbrlists)
!     ! SOLVE
!     WRITE(*,*) 'Solving heat conduction...'
!     CALL solve(blocks, nbrlists, iter, res_hist)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SAVE RESULTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    WRITE(*,*) 'Writing results...'

    !TURN THIS ON FOR PJ5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ! SAVE SOLUTION AS PLOT3D FILES
!     CALL plot3D(blocks)
    ! CALC TOTAL WALL TIME
    end_total = MPI_Wtime()
    wall_time_total = end_total - start_total

    !TURN THIS ON FOR PJ5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ! SAVE RESIDUAL HISTORY
!     CALL write_res(res_hist)
!     ! SAVE SOLVER PERFORMANCE PARAMETERS
!     CALL output(blocks, iter)

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

    WRITE(*,*) 'Done!'

    CALL MPI_Finalize(ierror)
END PROGRAM heatTrans
