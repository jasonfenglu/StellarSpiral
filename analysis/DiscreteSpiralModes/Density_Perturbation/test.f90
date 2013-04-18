program test
USE STELLARDISK_MODEL
USE STELLARDISK
DOUBLE PRECISION        ::ri,rf,dr,r
INTEGER                 ::n,i
type(typgalaxy_para)para1
type(typspiral)spiral
type tst
     integer,allocatable::dat(:)
endtype
type(tst) tt

ri = 0.d0
rf = 10.d0
!rf = 10.d0
N  = 1000

dr = (rf-ri)/dble(n)

!CALL stdpara.readstd
!CALL stdpara.printpara

!CALL spiral.init(spiral,100,12.d0,stdpara,1)
!CALL FindSpiral(spiral)
!CALL spiral.printh1
!DO i = 0, N
!        r = dr*dble(i)
!       write(6,*)r,sigma0(r,spiral)
!ENDDO


CALL stdpara.readstd
ALLOCATE(tt.dat(3))
!$OMP PARALLEL FIRSTPRIVATE(spiral)
!$OMP DO 
DO i = 1, 32
CALL spiral.init(spiral,500,12.d0,stdpara,1)
CALL FindSpiral(spiral)
CALL spiral.final
print *,spiral.error
ENDDO
!$OMP END PARALLEL

stop
end      
