program test
USE PLOTTING
USE STELLARDISK
USE OMP_LIB
USE STELLARDISK_MODEL
USE io
IMPLICIT NONE
include 'mpif.h'
type searchgrid_type
        DOUBLE PRECISION,ALLOCATABLE::coord(:,:,:)
        DOUBLE PRECISION,ALLOCATABLE::error(:,:)
        DOUBLE PRECISION,ALLOCATABLE::lcoord(:,:)
        DOUBLE PRECISION,ALLOCATABLE::lerror(:)
endtype
type(searchgrid_type)           ::searchgrid,recvgrid
type(typspiral)                 ::spiral
CHARACTER(len=32)                 ::arg
DOUBLE PRECISION                ::dr,wri,wii,di,err
DOUBLE PRECISION                ::domain(4) = (/40d0,100d0,0d0,-3d0/)  !better resolution 40-120,40 for 1 kpc
DOUBLE PRECISION                ::wr,wi
DOUBLE PRECISION,ALLOCATABLE    ::errormpisend(:),mpiall(:)
INTEGER                         ::l,i,j,p(1),n,m
INTEGER                         ::ipc
INTEGER                         ::now(3)
INTEGER                         ::complete_count = 0
INTEGER                         ::narg
INTEGER                         ::mpi_size,myid,info,ierr,chunk
INTEGER                         ::ompid, ompsize
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_size,ierr)
print *,'MPI started',myid


!reading hand put search range
narg = iargc()
SELECT CASE(narg)
CASE(2)
        CALL getarg(1,arg)
        READ(arg,*)domain(1)
        CALL getarg(2,arg)
        READ(arg,*)domain(2)
        IF(myid.eq.0)THEN
                print *,'only r defined:',domain
        ENDIF
CASE(4)
        DO i = 1, 4
                CALL getarg(i,arg)
                READ(arg,*)domain(i)
        ENDDO
        IF(myid.eq.0)THEN
                print *,'defined region:',domain
        ENDIF
CASE DEFAULT
ENDSELECT

IF(domain(2)-domain(1)<10.d0)THEN
        m = 400
ELSE
        m = int(domain(2)-domain(1))*40
ENDIF
n = m/10
if(mod(m,mpi_size).eq.0)then
        chunk = m/mpi_size
else 
        chunk = m/mpi_size + 1
endif
m = chunk*mpi_size
chunk = chunk *n

ALLOCATE(searchgrid.coord(m,n,2))
ALLOCATE(searchgrid.lcoord(m*n,2))

dr = (domain(2)-domain(1))/dble(m)
di = (domain(4)-domain(3))/dble(n)
!most left and upper grid
wri = domain(1)
wii = domain(3)
DO i = 1,m
        searchgrid%coord(i,:,1) = dble(i-1)*dr + wri
!       searchgrid%coord(i,:,1) =                wri
enddo
DO i = 1,n
        searchgrid%coord(:,i,2) = dble(i-1)*di + wii
!       searchgrid%coord(:,i,2) =                wii
enddo
searchgrid.lcoord = reshape(searchgrid.coord,(/m*n,2/))
DEALLOCATE(searchgrid.coord)

CALL stdpara.readstd
print *,myid,m,n,chunk
ALLOCATE(errormpisend(chunk))
!$OMP PARALLEL SHARED(searchgrid,complete_count,stdpara) FIRSTPRIVATE(spiral)
CALL spiral.init(100,15.d0,stdpara,1)
!$OMP DO 
!DO j = 1,m*n
DO j = chunk*myid+1,chunk*(myid+1)
        wr = searchgrid%lcoord(j,1)
        wi = searchgrid%lcoord(j,2)
        spiral.w = dcmplx(wr,wi)
        spiral.winit = .true.
!       ipc = omp_get_thread_num()
!       print *,'!!!',j,ipc
        CALL FindSpiral(spiral)
!       searchgrid%lerror(j) = abs(spiral.error)
        errormpisend(j-chunk*myid) = abs(spiral.error)
        if(myid.eq.0)then
                !$OMP CRITICAL
                complete_count = complete_count + 1
                print *,real(complete_count)/real(chunk)*100.
                !$OMP END CRITICAL
        endif
ENDDO
!$OMP END DO
CALL spiral.free
!$OMP END PARALLEL
print *,'collecting data',myid

ALLOCATE(searchgrid.lerror(m*n))
CALL MPI_GATHER(errormpisend,chunk,MPI_DOUBLE_PRECISION,searchgrid.lerror,chunk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
DEALLOCATE(errormpisend)
ALLOCATE(searchgrid.error(m,n))

!p = MINLOC(searchgrid%lerror(:))
!wri = searchgrid%lcoord(p(1),1)
!wii = searchgrid%lcoord(p(1),2)
!err = searchgrid%lerror(p(1))
!DO i = 1, N
!DO j = 1, N
!        print *,searchgrid.coord(i,j,:),searchgrid.error(i,j)
!ENDDO
!ENDDO

if(myid.eq.0)then
        searchgrid.error = reshape(searchgrid.lerror,(/m,n/))
        print *,'min error',minval(searchgrid.error(:,:))
        CALL plotpspdsearch(searchgrid.error,m,n,domain)
endif

!1000 CALL INIT_STELLARDISK(100,40.d0)
!wr = 60.d0
!wi = -1.d0
!CALL FindSpiral
!DO i = 1, 100
!        write(10,*)spiral.r(i),real(k3sqrt(spiral.r(i)))
!enddo
!DEALLOCATE(searchgrid.coord)
!DEALLOCATE(searchgrid.error)
!DEALLOCATE(searchgrid.lcoord)
!DEALLOCATE(searchgrid.lerror)

CALL MPI_FINALIZE(ierr)
ENDPROGRAM
