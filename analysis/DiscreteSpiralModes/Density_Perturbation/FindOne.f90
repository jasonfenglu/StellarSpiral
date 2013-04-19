PROGRAM find_all
USE STELLARDISK_MODEL
IMPLICIT NONE
type(typspiral),TARGET            ::spiral0
DOUBLE PRECISION                  ::wri,wii
CHARACTER(len=32)                 ::arg
INTEGER                           ::mode,narg

narg = iargc()
SELECT case(narg)
CASE(2)
        CALL getarg(1,arg)
        READ(arg,*)wri
        CALL getarg(2,arg)
        READ(arg,*)wii
CASE(1)
        print *,'use input mode'
        CALL getarg(1,arg)
        READ(arg,*)mode
!       CALL spiral0.init(spiral0,100,12.d0,stdpara,mode)
        CALL spiral0.readw(mode)
        wri = real(spiral0.w)
        wii = imag(spiral0.w)
!       CALL spiral0.final
CASE DEFAULT        
        print *,'no input initial finding value, using default.'
!       CALL spiral0.init(spiral0,100,12.d0,stdpara,1)
        CALL spiral0.readw(1)
        wri = real(spiral0.w)
        wii = imag(spiral0.w)
!       CALL spiral0.final
ENDSELECT

          
CALL findpspsd(wri,wii)
print *,wri,',',wii

STOP
END PROGRAM

SUBROUTINE omp_single_grid(l,wri,wii,err,r)
USE STELLARDISK_MODEL
USE STELLARDISK,only:FindSpiral,pi_n=>pi
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
type(typspiral)                 ::spiral
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

CALL stdpara.readstd
!CALL omp_set_num_threads(1)
!$OMP PARALLEL SHARED(searchgrid,stdpara) FIRSTPRIVATE(spiral)
!$OMP DO 
DO j = 1,n**2
        CALL spiral.init(spiral,200,12.d0,stdpara,1)
        spiral.w = dcmplx(searchgrid.lcoord(j,1),searchgrid.lcoord(j,2))
        CALL FindSpiral(spiral)
        searchgrid.lerror(j) = abs(spiral.error)
        CALL spiral.final
ENDDO
!$OMP END DO 
!$OMP END PARALLEL 
searchgrid.error = reshape(searchgrid.lerror,(/n,n/))
p = MINLOC(searchgrid%error(:,:))
i = p(1)
j = p(2)
wri = searchgrid%coord(i,j,1)
wii = searchgrid%coord(i,j,2)
err = searchgrid%error(i,j)

if(j.eq.1 .or. j.eq.n .or. i.eq.1 .or. i.eq.n)then
        l = l - 2
        CALL XERMSG('k3sqrt','Omega Finding','Eigenvalue at boundary, using &
        larger grid.',-95,-1)
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
IMPLICIT NONE
DOUBLE PRECISION                ::wri,wii,err,r
INTEGER                         ::l
l = 1
write(*,'(I2,3X,F7.4,3X,F7.4,3X,E10.3)')0,wri,wii
do while (l.le.20)
        CALL omp_single_grid(l,wri,wii,err,r)
        write(*,'(I2,3X,F7.4,3X,F7.4,3X,E10.3,3X,D10.3,3X,F7.4)')l,wri,wii,err,r
        if(wii.gt.0.d0)then
                write(0,*)"don't growth, stop"
                stop
        endif
        if(abs(err).le.1d-6)exit
        l = l + 1
enddo

ENDSUBROUTINE 
