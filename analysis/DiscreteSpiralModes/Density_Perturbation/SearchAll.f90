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
DOUBLE PRECISION                ::domain(4) = (/30d0,60d0,0d0,-4d0/)
INTEGER                         ::l,i,j,p(1),n,m
INTEGER                         ::ipc
INTEGER                         ::now(3)
INTEGER                         ::complete_count = 0

m = 800
n = m/5

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

!$OMP PARALLEL SHARED(searchgrid,complete_count) PRIVATE(spiral,stdpara)
CALL INIT_STELLARDISK(200,15.d0)
!$OMP DO PRIVATE(spiral)
DO j = 1,m*n
        wr = searchgrid%lcoord(j,1)
        wi = searchgrid%lcoord(j,2)
!       ipc = omp_get_thread_num()
!       print *,'!!!',j,ipc
        CALL Findu
        searchgrid%lerror(j) = abs(error())
        !$OMP CRITICAL
                complete_count = complete_count + 1
        !$OMP END CRITICAL
        print *,real(complete_count)/real(m*n)*100.
ENDDO
!$OMP END DO
CALL ENDSTELLARDISK
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
