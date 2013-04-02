PROGRAM find_all
USE STELLARDISK,ONLY:INIT_STELLARDISK,ENDSTELLARDISK,wr,wi
IMPLICIT NONE
DOUBLE PRECISION                  ::wri,wii
CHARACTER(len=32)                 ::arg
          
if(iargc().eq.2)then
        CALL getarg(1,arg)
        READ(arg,*)wri
        CALL getarg(2,arg)
        READ(arg,*)wii
else 
        CALL INIT_STELLARDISK(5,20.d0)
        print *,'no input initial finding value, using default'
        wri = wr
        wii = wi
        CALL ENDSTELLARDISK
endif
CALL findpspsd(wri,wii)
print *,wri,',',wii

STOP
END PROGRAM

SUBROUTINE single_grid(l,wri,wii)
USE STELLARDISK
IMPLICIT NONE
type searchgrid_type
sequence
        DOUBLE PRECISION::coord(12,12,2)
        DOUBLE PRECISION::error(12,12)
endtype
type(searchgrid_type)           ::searchgrid,recvgrid
DOUBLE PRECISION                ::dr,wri,wii,di
INTEGER                         ::l,i,j,p(2)



dr = 1.d0/10.0d0**(l-1)
di = 0.5d0/10.0d0**(l-1)
!most left and upper grid
wri = wri +(-6.d0+0.5d0)*dr
wii = wii +(-6.d0+0.5d0)*dr

DO i = 1,12
        searchgrid%coord(:,i,2) = dble(i-1)*dr + wii
        searchgrid%coord(i,:,1) = dble(i-1)*dr + wri
enddo
DO j = 1,12
CALL INIT_STELLARDISK(500,20.d0)
DO i = 1,12
        wr = searchgrid%coord(i,j,1)
        wi = searchgrid%coord(i,j,2)
        CALL FindSpiral
ENDDO
ENDDO
        
p = MINLOC(searchgrid%error(:,:))
i = p(1)
j = p(2)
wri = searchgrid%coord(i,j,1)
wii = searchgrid%coord(i,j,2)

!!Print search grid for debug
!if(j.eq.1 .or. j.eq.12 .or. i.eq.1 .or. i.eq.12)l = l -1
!DO i = 1,12
!DO j = 1,12
!        print *,searchgrid%coord(i,j,:),searchgrid%error(i,j)
!ENDDO
!ENDDO

ENDSUBROUTINE

SUBROUTINE omp_single_grid(l,wri,wii,err,r)
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
DOUBLE PRECISION                ::r
INTEGER                         ::n
INTEGER                         ::l,i,j,p(2)
INTEGER                         ::ipc
INTEGER                         ::now(3)

r  = 2.d0**(1.d0-dble(l))
n  = 8
dr = r/dble(n)
di = dr
ALLOCATE(searchgrid.coord(n,n,2))
ALLOCATE(searchgrid.error(n,n))
ALLOCATE(searchgrid.lcoord(n*n,2))
ALLOCATE(searchgrid.lerror(n*n))

!most left and upper grid
wri  = wri - r/2.d0
wii  = wii - r/2.d0
!wri = wri +(-r/2.d0+0.5d0)*dr
!wii = wii +(-r/2.d0+0.5d0)*di
DO i = 1,n
        searchgrid%coord(:,i,2) = dble(i-1)*dr + wii
        searchgrid%coord(i,:,1) = dble(i-1)*dr + wri
!       searchgrid%coord(:,i,2) =                wii
!       searchgrid%coord(i,:,1) =                wri
enddo

searchgrid.lcoord = reshape(searchgrid.coord,(/n*n,2/))

!$OMP PARALLEL SHARED(searchgrid) PRIVATE(spiral,stdpara)
CALL INIT_STELLARDISK(200,13.d0)
!$OMP DO PRIVATE(spiral)
DO j = 1,n**2
        wr = searchgrid%lcoord(j,1)
        wi = searchgrid%lcoord(j,2)
!       ipc = omp_get_thread_num()
!       print *,'!!!',j,ipc
        CALL Findu
        searchgrid%lerror(j) = abs(error())
ENDDO
!$OMP END DO
!$OMP END PARALLEL
CALL ENDSTELLARDISK
searchgrid.error = reshape(searchgrid.lerror,(/n,n/))
p = MINLOC(searchgrid%error(:,:))
i = p(1)
j = p(2)
wri = searchgrid%coord(i,j,1)
wii = searchgrid%coord(i,j,2)
err = searchgrid%error(i,j)

if(j.eq.1 .or. j.eq.n .or. i.eq.1 .or. i.eq.n)then
        l = l - 2
        CALL XERMSG('k3sqrt','Omega Finding','Eigenvalue at boundary, giving up.',-95,-1)
endif
!DO i = 1,12
!DO j = 1,12
!        print *,searchgrid%coord(i,j,:),searchgrid%error(i,j)
!ENDDO
!ENDDO
DEALLOCATE(searchgrid.coord)
DEALLOCATE(searchgrid.error)
DEALLOCATE(searchgrid.lcoord)
DEALLOCATE(searchgrid.lerror)
ENDSUBROUTINE

SUBROUTINE findpspsd(wri,wii)
USE STELLARDISK,only:spiral
IMPLICIT NONE
DOUBLE PRECISION                ::wri,wii,err,r
INTEGER                         ::l
l = 1
write(*,'(I2,3X,F7.4,3X,F7.4,3X,E10.3)')0,wri,wii
do while (l.le.20)
        CALL omp_single_grid(l,wri,wii,err,r)
        write(*,'(I2,3X,F7.4,3X,F7.4,3X,E10.3,3X,D10.3,3X,F7.4)')l,wri,wii,err,r,spiral.fortoone
        if(abs(err).le.1d-6)exit
        l = l + 1
enddo

ENDSUBROUTINE 
