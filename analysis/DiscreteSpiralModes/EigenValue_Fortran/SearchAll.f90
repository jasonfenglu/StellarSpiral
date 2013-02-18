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
DOUBLE PRECISION                ::domain(4) = (/15d0,65d0,0d0,-2d0/)
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
!searchgrid.lerror = reshape(searchgrid.error,(/144/))

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
DO i = 1, N
DO j = 1, N
        print *,searchgrid.coord(i,j,:),searchgrid.error(i,j)
ENDDO
ENDDO



!if(j.eq.1 .or. j.eq.12 .or. i.eq.1 .or. i.eq.12)l = l -1
!DO i = 1,12
!DO j = 1,12
!        print *,searchgrid%coord(i,j,:),searchgrid%error(i,j)
!ENDDO
!ENDDO

CALL plot2d(searchgrid.error,n,n,domain)

100 CALL INIT_STELLARDISK(100,40.d0)
wr = 10
wi = 10
CALL FindSpiral
DO i = 1, 100
        write(10,*)spiral.r(i),real(Sigma0(spiral.r(i)))
enddo


endprogram
