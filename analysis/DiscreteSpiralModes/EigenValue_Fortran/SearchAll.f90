program test
USE PLOTTING
USE STELLARDISK
USE OMP_LIB

IMPLICIT NONE
type searchgrid_type
sequence
        DOUBLE PRECISION,ALLOCATABLE::coord(:,:,:)
        DOUBLE PRECISION,ALLOCATABLE::error(:,:)
        DOUBLE PRECISION,ALLOCATABLE::lcoord(:,:)
        DOUBLE PRECISION,ALLOCATABLE::lerror(:)
endtype
type(searchgrid_type)           ::searchgrid,recvgrid
DOUBLE PRECISION                ::dr,wri,wii,di,err
DOUBLE PRECISION                ::domain(4) = (/15d0,70d0,0d0,-4d0/)
INTEGER                         ::l,i,j,p(1),n
INTEGER                         ::ipc
INTEGER                         ::now(3)

n = 100
ALLOCATE(searchgrid.coord(n,n,2))
ALLOCATE(searchgrid.error(n,n))
ALLOCATE(searchgrid.lcoord(n*n,2))
ALLOCATE(searchgrid.lerror(n*n))

dr = (domain(2)-domain(1))/dble(n)
di = (domain(4)-domain(3))/dble(n)
!most left and upper grid
wri = domain(1)
wii = domain(3)
DO i = 1,n
        searchgrid%coord(i,:,1) = dble(i-1)*dr + wri
        searchgrid%coord(:,i,2) = dble(i-1)*di + wii
!       searchgrid%coord(i,:,1) =                wri
!       searchgrid%coord(:,i,2) =                wii
enddo
searchgrid.lcoord = reshape(searchgrid.coord,(/n*n,2/))
!$OMP PARALLEL 
CALL INIT_STELLARDISK(200,20.d0)
!$OMP BARRIER
!$OMP DO ORDERED
DO j = 1,n*n
        wr = searchgrid%lcoord(j,1)
        wi = searchgrid%lcoord(j,2)
!       ipc = omp_get_thread_num()
!       print *,'!!!',j,ipc
        CALL FindSpiral
        searchgrid%lerror(j) = abs(error())
ENDDO
!$OMP END DO
!$OMP END PARALLEL
CALL ENDSTELLARDISK
!p = MINLOC(searchgrid%lerror(:))
!wri = searchgrid%lcoord(p(1),1)
!wii = searchgrid%lcoord(p(1),2)
!err = searchgrid%lerror(p(1))
searchgrid.coord = reshape(searchgrid.coord,(/n,n,2/))
searchgrid.error = reshape(searchgrid.lerror,(/n,n/))
!DO i = 1, N
!DO j = 1, N
!        print *,searchgrid.coord(i,j,:),searchgrid.error(i,j)
!ENDDO
!ENDDO


CALL plot2d(searchgrid.error,n,n,domain)

1000 CALL INIT_STELLARDISK(100,40.d0)
wr = 60.d0
wi = -1.d0
CALL FindSpiral
DO i = 1, 100
        write(10,*)spiral.r(i),real(k3sqrt(spiral.r(i)))
enddo


endprogram
