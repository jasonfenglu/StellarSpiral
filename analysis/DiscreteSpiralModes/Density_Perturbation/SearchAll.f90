program test
USE PLOTTING
USE STELLARDISK
USE OMP_LIB
USE STELLARDISK_MODEL

IMPLICIT NONE
type searchgrid_type
sequence
        DOUBLE PRECISION,ALLOCATABLE::coord(:,:,:)
        DOUBLE PRECISION,ALLOCATABLE::error(:,:)
        DOUBLE PRECISION,ALLOCATABLE::lcoord(:,:)
        DOUBLE PRECISION,ALLOCATABLE::lerror(:)
endtype
type(searchgrid_type)           ::searchgrid,recvgrid
type(typspiral)                 ::spiral
CHARACTER(len=32)                 ::arg
DOUBLE PRECISION                ::dr,wri,wii,di,err
DOUBLE PRECISION                ::domain(4) = (/40d0,120d0,0d0,-3d0/)  !better resolution 40-120
DOUBLE PRECISION                ::wr,wi
INTEGER                         ::l,i,j,p(1),n,m
INTEGER                         ::ipc
INTEGER                         ::now(3)
INTEGER                         ::complete_count = 0
INTEGER                         ::narg

narg = iargc()

SELECT CASE(narg)
CASE(2)
                CALL getarg(1,arg)
                READ(arg,*)domain(1)
                CALL getarg(2,arg)
                READ(arg,*)domain(2)
CASE DEFAULT
ENDSELECT

!m = 800
m = int(domain(2)-domain(1))*20
n = m/10

ALLOCATE(searchgrid.coord(m,n,2))
ALLOCATE(searchgrid.error(m,n))
ALLOCATE(searchgrid.lcoord(m*n,2))
ALLOCATE(searchgrid.lerror(m*n))

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

CALL stdpara.readstd
!$OMP PARALLEL SHARED(searchgrid,complete_count,stdpara) FIRSTPRIVATE(spiral)
!$OMP DO 
DO j = 1,m*n
        CALL spiral.init(spiral,100,15.d0,stdpara,1)
        wr = searchgrid%lcoord(j,1)
        wi = searchgrid%lcoord(j,2)
        spiral.w = dcmplx(wr,wi)
!       ipc = omp_get_thread_num()
!       print *,'!!!',j,ipc
        CALL FindSpiral(spiral)
        searchgrid%lerror(j) = abs(spiral.error)
        CALL spiral.final
        !$OMP CRITICAL
                complete_count = complete_count + 1
        !$OMP END CRITICAL
        print *,real(complete_count)/real(m*n)*100.
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!p = MINLOC(searchgrid%lerror(:))
!wri = searchgrid%lcoord(p(1),1)
!wii = searchgrid%lcoord(p(1),2)
!err = searchgrid%lerror(p(1))
searchgrid.coord = reshape(searchgrid.coord,(/m,n,2/))
searchgrid.error = reshape(searchgrid.lerror,(/m,n/))
!DO i = 1, N
!DO j = 1, N
!        print *,searchgrid.coord(i,j,:),searchgrid.error(i,j)
!ENDDO
!ENDDO
print *,'min error',minval(searchgrid.error(:,:))


CALL plotpspdsearch(searchgrid.error,m,n,domain)

!1000 CALL INIT_STELLARDISK(100,40.d0)
!wr = 60.d0
!wi = -1.d0
!CALL FindSpiral
!DO i = 1, 100
!        write(10,*)spiral.r(i),real(k3sqrt(spiral.r(i)))
!enddo
DEALLOCATE(searchgrid.coord)
DEALLOCATE(searchgrid.error)
DEALLOCATE(searchgrid.lcoord)
DEALLOCATE(searchgrid.lerror)

ENDPROGRAM
