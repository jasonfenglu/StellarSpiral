MODULE PspdSearchStruct
type searchgrid_type
sequence
        DOUBLE PRECISION,ALLOCATABLE::coord(:,:,:)
        DOUBLE PRECISION,ALLOCATABLE::error(:,:)
        DOUBLE PRECISION,ALLOCATABLE::lcoord(:,:)
        DOUBLE PRECISION,ALLOCATABLE::lerror(:)
endtype

ENDMODULE

program test
USE PspdSearchStruct
USE STELLARDISK
type(searchgrid_type),ALLOCATABLE  ::searchset(:)
DOUBLE PRECISION,TARGET,ALLOCATABLE::paraset(:,:)
DOUBLE PRECISION,ALLOCATABLE       ::packtoplot(:,:,:)
DOUBLE PRECISION                   ::domain(4) = (/40d0,200,0d0,-6d0/)
INTEGER                            ::m,n
INTEGER                            ::nmodel = 9
INTEGER                            ::i,j,k


m = 2000
n = m/5

!Read in Standard Parameters
CALL READINSTDPARA
ALLOCATE(paraset(size(stdpara,1),nmodel))
DO i = 1, nmodel
        paraset(:,i) = stdpara(:)
ENDDO
!Setting deviations
DO i = -4,4,1
        paraset(10,i+5) = paraset(10,5)*(1.d0+dble(i)*1d-1)
ENDDO

!Distributed Calculation
ALLOCATE(searchset(nmodel))
DO i = 1, nmodel
        CALL SearchOneModel(searchset(i),i)
ENDDO
!Collect data to plot subroutine
ALLOCATE(packtoplot(m,n,nmodel))
DO i = 1, nmodel
        packtoplot(:,:,i) = searchset(i).error
ENDDO
CALL plotpspdsearchset(packtoplot,m,n,domain)

CONTAINS
SUBROUTINE SearchOneModel(resultgrid,k)
USE PLOTTING
USE STELLARDISK
USE OMP_LIB
USE PspdSearchStruct
IMPLICIT NONE

type(searchgrid_type)           ::searchgrid,resultgrid
DOUBLE PRECISION                ::dr,wri,wii,di,err
INTEGER                         ::l,i,j,p(1),k
INTEGER                         ::ipc
INTEGER                         ::now(3)
INTEGER                         ::complete_count = 0


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

!$OMP PARALLEL SHARED(searchgrid,complete_count,para) PRIVATE(spiral)
CALL INIT_STELLARDISK(200,15.d0)
para=>paraset(:,k)
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


!CALL plotpspdsearch(searchgrid.error,m,n,domain)

!1000 CALL INIT_STELLARDISK(100,40.d0)
!wr = 60.d0
!wi = -1.d0
!CALL FindSpiral
!DO i = 1, 100
!        write(10,*)spiral.r(i),real(k3sqrt(spiral.r(i)))
!enddo
resultgrid = searchgrid
DEALLOCATE(searchgrid.coord)
DEALLOCATE(searchgrid.error)
DEALLOCATE(searchgrid.lcoord)
DEALLOCATE(searchgrid.lerror)


ENDSUBROUTINE

ENDPROGRAM
