! MAE 267
! PROJECT 5
! LOGAN HALSTROM
! 29 NOVEMBER 2015

! DESCRIPTION:  Subroutines used for solving heat conduction of steel plate.
! Subroutines utilizing linked lists are here so that linked lists do not need
! to be function inputs.
! Utilizes modules from 'modules.f90'

! CONTENTS:
    ! init_gridsystem
        ! Initialize the solution with dirichlet B.C.s.  Save to restart files.

    ! init_solution
        ! Read initial conditions from restart files.  Then calculate parameters
        ! used in solution

    ! solve
        ! Solve heat conduction equation with finite volume scheme
        ! (within iteration loop)

    ! output
        ! Save solution performance parameters to file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE subroutines
    USE CONSTANTS
    USE BLOCKMOD
    USE IO

    IMPLICIT NONE

    ! SOLUTION BLOCKS
    ! (initialized individually for each parallel processor,
    !   holds specific blocks distributed to each specific processor)
    TYPE(BLKTYPE), POINTER :: blocks(:)
    ! LINKED LISTS STORING NEIGHBOR INFO
    TYPE(NBRLIST) :: nbrlists
    ! neighbors on other processors
    TYPE(NBRLIST) :: mpilists

CONTAINS
    SUBROUTINE init_gridsystem(blocks, procs)
        ! Initialize the solution with dirichlet B.C.s.  Save to restart files.

        TYPE(BLKTYPE)  :: blocks(:)
        TYPE(PROCTYPE) :: procs(:)

        ! INITIALIZE BLOCKS
        CALL init_blocks(blocks)

        ! CALC LOCAL BOUNDARIES OF CELLS
        CALL set_block_bounds(blocks)

        ! INITIALIZE MESH
        CALL init_mesh(blocks)
        ! INITIALIZE TEMPERATURE WITH DIRICHLET B.C.
        CALL init_temp(blocks)

        ! DISTRIBUTE BLOCKS TO PROCESSORS
        CALL dist_blocks(blocks, procs)
        ! DETERMIN NEIGHBOR PROCESSOR INFORMATION
        CALL init_neighbor_procs(blocks, procs)


        ! WRITE BLOCK CONNECTIVITY FILE
        CALL write_config(procs)

    END SUBROUTINE init_gridsystem

    SUBROUTINE init_solution(blocks, nbrlists, mpilists)
        ! Read initial conditions from restart files.  Then calculate parameters
        ! used in solution

        TYPE(BLKTYPE), POINTER  :: blocks(:)
        ! LINKED LISTS STORING NEIGHBOR INFO
        TYPE(NBRLIST) :: nbrlists, mpilists

        write(*,*) "read config", MYID
        ! READ BLOCK CONFIGURATION INFORMATION FROM CONFIG FILE
        CALL read_config(blocks)

        ! INITIALIZE LINKED LISTS CONTAINING BOUNDARY INFORMATION
        write(*,*) 'make linked lists', MYID
        CALL init_linklists(blocks, nbrlists, mpilists)
        ! POPULATE BLOCK GHOST NODES
        write(*,*) 'update ghosts', MYID
        CALL update_ghosts_sameproc(blocks, nbrlists)
        CALL update_ghosts_diffproc_send(blocks, mpilists)
        CALL update_ghosts_diffproc_recv(blocks, mpilists)

        ! CALC AREAS FOR SECONDARY FLUXES
        write(*,*) 'calc solution stuff', MYID
        CALL calc_cell_params(blocks)
        ! CALC CONSTANTS OF INTEGRATION
        CALL calc_constants(blocks)
        write(*,*) 'calced solution stuff', MYID

    END SUBROUTINE init_solution


    SUBROUTINE solve(blocks, nbrlists, mpilists, iter, res_hist)
        ! Solve heat conduction equation with finite volume scheme
        ! (within iteration loop)

        TYPE(BLKTYPE) :: blocks(:)
        ! LINKED LISTS STORING NEIGHBOR INFO
        TYPE(NBRLIST) :: nbrlists, mpilists
        ! Residual history linked list
        TYPE(RESLIST), POINTER :: res_hist
        ! pointer to iterate linked list
        TYPE(RESLIST), POINTER :: hist
        ! Minimum residual criteria for iteration, actual residual
        REAL(KIND=8) :: res = 1000.D0, resloc=0.D0, resmax=0.D0
        ! iter in function inputs so it can be returned to main
        INTEGER :: iter, IBLK, IBLKRES


        REAL(KIND=8) :: start_solve, end_solve
        IF (MYID == 0) THEN
            ! START SOLVER CLOCK
            start_solve = MPI_Wtime()
        END IF

        ! residual history
        ALLOCATE(res_hist)
        hist => res_hist

        iter_loop: DO WHILE (res >= min_res .AND. iter <= max_iter)
            ! Iterate FV solver until residual becomes less than cutoff or
            ! iteration count reaches given maximum

            ! CALC NEW TEMPERATURE AT ALL POINTS
            CALL calc_temp(blocks)

            ! UPDATE GHOST NODES WITH NEW TEMPERATURE SOLUTION
            CALL update_ghosts_sameproc(blocks, nbrlists)
            CALL update_ghosts_diffproc_send(blocks, mpilists)
            CALL update_ghosts_diffproc_recv(blocks, mpilists)

            ! CALC RESIDUAL FOR LOCAL BLOCKS
            resmax = 0.D0
            DO IBLK = 1, MYNBLK
                ! Find max of each block
                resloc = MAXVAL( ABS( blocks(IBLK)%mesh%Ttmp(2:IMAXBLK-1, 2:JMAXBLK-1) ) )
                ! keep biggest residual
                IF (resloc > resmax) THEN
                    resmax = resloc
                END IF
            END DO
!             write(*,*) "before mpiallreduce:", MYID, resmax
            ! FINAL MAX RESIDUAL (FOR ALL PROCESSORS)
            CALL MPI_ALLREDUCE(resmax, res, 1, MPI_REAL8, MPI_MAX, &
                                    MPI_COMM_WORLD, IERROR)
!             write(*,*) "after mpiallreduce:", MYID, res

            ! SWITCH TO NEXT LINK
                ! (skip first entry)
            ALLOCATE(hist%next)
            hist => hist%next
            NULLIFY(hist%next)
            ! STORE RESIDUAL HISTORY
            hist%iter = iter
            hist%res = res


            ! INCREMENT ITERATION COUNT
            iter = iter + 1

        END DO iter_loop

        ! HOLD UNTIL ALL PROCESSORS HAVE FINISHED ITERATION LOOP
        CALL MPI_Barrier(MPI_COMM_WORLD, IERROR)

        ! there was an extra increment after final iteration we need to subtract
        iter = iter - 1

        IF (MYID == 0) THEN
            ! CALC SOLVER WALL CLOCK TIME
            end_solve = MPI_Wtime()
            wall_time_solve = end_solve - start_solve
        END IF

        IF (iter > max_iter) THEN
          WRITE(*,*) 'DID NOT CONVERGE (NUMBER OF ITERATIONS:', iter, ')'
        ELSE
          WRITE(*,*) 'CONVERGED (NUMBER OF ITERATIONS:', iter, ')'
          WRITE(*,*) '          (MAXIMUM RESIDUAL    :', res,  ')'
        END IF
    END SUBROUTINE solve

    SUBROUTINE output(blocks, iter)
        ! Save solution performance parameters to file

        TYPE(BLKTYPE), TARGET :: blocks(:)
        REAL(KIND=8), POINTER :: tmpT(:,:), tempTemperature(:,:)
        REAL(KIND=8) :: resloc, resmax
        INTEGER :: iter, I, J, IBLK, IRES

!         Temperature => mesh%T(2:IMAX-1, 2:JMAX-1)
!         tempTemperature => mesh%Ttmp(2:IMAX-1, 2:JMAX-1)

        ! CALC RESIDUAL
        resmax = 0.D0
        DO IBLK = 1, MYNBLK
            ! Find max of each block
            resloc = MAXVAL( ABS( blocks(IBLK)%mesh%Ttmp(2:IMAXBLK-1, 2:JMAXBLK-1) ) )
            ! keep biggest residual
            IF (resloc > resmax) THEN
                resmax = resloc
                IRES = IBLK
            END IF
        END DO

        CALL MPI_Barrier(MPI_COMM_WORLD, IERROR)
!         CALL MPI_Bcast(MPI_COMM_WORLD, IERROR)


        ! Write final maximum residual and location of max residual
!         OPEN(UNIT = 1, FILE = casedir // "SteadySoln.dat")
!         DO i = 1, IMAX
!             DO j = 1, JMAX
!                 WRITE(1,'(F10.7, 5X, F10.7, 5X, F10.7, I5, F10.7)'), mesh%x(i,j), mesh%y(i,j), mesh%T(i,j)
!             END DO
!         END DO
!         CLOSE (1)

        ! Screen output
        tmpT => blocks(IRES)%mesh%Ttmp
        WRITE (*,*), "IMAX/JMAX", IMAX, JMAX
        WRITE (*,*), "N/M", N, M
        WRITE (*,*), "iters", iter
        WRITE (*,*), "max residual", MAXVAL(tmpT(2:IMAXBLK-1, 2:JMAXBLK-1))
        WRITE (*,*), "on block id", IRES
        WRITE (*,*), "residual ij", MAXLOC(tmpT(2:IMAXBLK-1, 2:JMAXBLK-1))

        ! Write to file
        IF (MYID == 0) THEN
    !         OPEN (UNIT = 2, FILE = TRIM(casedir) // "SolnInfo.dat")
            OPEN (UNIT = 2, FILE = "SolnInfo.dat")
            WRITE (2,*), "Running a", IMAX, "by", JMAX, "grid,"
            WRITE (2,*), "With NxM:", N, "x", M, "blocks took:"
            WRITE (2,*), iter, "iterations"
!             WRITE (2,*), wall_time_total, "seconds (Total CPU walltime)"
            WRITE (2,*), wall_time_solve, "seconds (Solver CPU walltime)"
    !         WRITE (2,*), wall_time_iter, "seconds (Iteration CPU walltime)"
            WRITE (2,*)
            WRITE (2,*), "Found max residual of ", MAXVAL(tmpT(2:IMAXBLK-1, 2:JMAXBLK-1))
            WRITE (2,*), "on block id", IRES
            WRITE (2,*), "At ij of ", MAXLOC(tmpT(2:IMAXBLK-1, 2:JMAXBLK-1))
            CLOSE (2)
        END IF
    END SUBROUTINE output



END MODULE subroutines


