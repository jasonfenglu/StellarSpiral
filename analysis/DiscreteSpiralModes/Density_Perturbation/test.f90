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
rf = 15.d0
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


!CALL stdpara.readstd
!ALLOCATE(tt.dat(3))
!!$OMP PARALLEL FIRSTPRIVATE(spiral)
!!$OMP DO 
!DO i = 1, 32
!CALL spiral.init(spiral,500,12.d0,stdpara,1)
!CALL FindSpiral(spiral)
!CALL spiral.final
!print *,spiral.error
!ENDDO
!!$OMP END PARALLEL

CALL stdpara.readstd
CALL spiral.init(spiral,N,12.d0,stdpara,1)
CALL spiral.readw(1)
CALL FindSpiral(spiral)
DO i = 1, N
        r = dr*dble(i)
!      write(6,*)r,real(spiral.u(2,i)),imag(spiral.u(2,i))
!      write(6,*)r,real(k3sqrt(r,spiral)),imag(k3sqrt(r,spiral))
        write(6,*)r,F2(r)
        
ENDDO

stop
CONTAINS

FUNCTION F(r)
IMPLICIT NONE
DOUBLE PRECISION                ::r,F
F = sqrt(k3sqrt(r,spiral))
ENDFUNCTION
Function F2(r) result(ans)
USE NUM
        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(IN)     ::r
        DOUBLE PRECISION                ::a,rr
        DOUBLE PRECISION                ::ans
        DOUBLE PRECISION                ::ferr = 1.d-15
        INTEGER                         ::IERR,K=6000
        a = zerolimit
        rr = r
        IERR = 0
        if(rr.eq.0.d0)then
                ans = 0.d0
        else
                CALL DGAUS8(F,a,rr,fERR,ans,IERR)
!               CALL DQNC79(F,a,rr,fERR,ans,IERR,K)
        endif

ENDFUNCTION 

end      
