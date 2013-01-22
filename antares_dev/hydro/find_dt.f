        MODULE find_dt

        type :: dt_struct
        sequence
                INTEGER         ::id
                INTEGER         ::direction
                INTEGER         ::i
                INTEGER         ::j
                REAL*8          ::dt
        endtype dt_struct

        type(dt_struct), SAVE   ::dts
        CONTAINS

        SUBROUTINE find_dt_init()
        dts = dt_struct(0,0,0,0,0.d0)
        ENDSUBROUTINE find_dt_init

        SUBROUTINE exchange_dts()
        USE COMMON_PARAMS,ONLY: nprocs,myid
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        include 'mpif.h'
        type(dt_struct),ALLOCATABLE     :: dtss(:)
        integer dt_type, oldtypes(0:1), blockcounts(0:1), 
     &        offsets(0:1), extent
        REAL*8,ALLOCATABLE              :: dtts(:)
        INTEGER                         ::barrier_com

        ALLOCATE(dtss(0:nprocs-1))
        ALLOCATE(dtts(0:nprocs-1))

ccccc  construct mpi package for dts
        !id, direction, i, j
        offsets(0) = 0
        oldtypes(0) = MPI_INTEGER
        blockcounts(0) = 4
        !dt
        call MPI_TYPE_EXTENT(MPI_REAL, extent, ierr)
        offsets(1) = 4 * extent
        oldtypes(1) = MPI_DOUBLE_PRECISION
        blockcounts(1) = 1
        !commit type for mpi
        call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, 
     &                       dt_type, ierr)
        call MPI_TYPE_COMMIT(dt_type, ierr)

ccccc  copy all dts to id=0
        !send mines
        if (myid.ne.0) then
              call MPI_Send(dts,1,dt_type,0,1,
     &        MPI_COMM_WORLD,ierr)
        endif
        !get others
        if (myid.eq.0) then
              dtss(0) = dts

              !recive
              do I=1,nprocs-1
              call MPI_Recv(dtss(I),1,dt_type,I,1,
     &        MPI_COMM_WORLD,istatus,ierr)
              enddo

              !searching min dt
              dts = dtss(MINLOC(dtss(:)%dt,1)-1)


        endif
        !send dts to everyone and to recieve new dts
        call MPI_Bcast(dts,1,dt_type,0,MPI_COMM_WORLD,ierr)

        call MPI_TYPE_FREE(dt_type, ierr)
        DEALLOCATE(dtss)
        DEALLOCATE(dtts)
c       write(6,*) 'LLSC complete myid=',myid
        ENDSUBROUTINE exchange_dts

        ENDMODULE find_dt
