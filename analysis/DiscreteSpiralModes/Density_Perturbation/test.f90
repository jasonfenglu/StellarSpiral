program test
USE STELLARDISK_MODEL
USE STELLARDISK
USE RK
DOUBLE PRECISION        ::ri,rf,dr,r
DOUBLE COMPLEX          ::h1
INTEGER                 ::n,i
type(typgalaxy_para)para1
type(typspiral)spiral
type tst
     integer,allocatable::dat(:)
endtype
type(tst) tt
DOUBLE COMPLEX          ::Y(2,20),Y1(2)

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

!CALL stdpara.readstd
!CALL spiral.init(N,12.d0,stdpara,1)
!CALL spiral.readw(1)
!CALL FindSpiral(spiral)
!DO i = 1, N
!        r = dr*dble(i)
!!      write(6,*)r,real(spiral.u(2,i)),imag(spiral.u(2,i))
!!      write(6,*)r,real(k3sqrt(r,spiral)),imag(k3sqrt(r,spiral))
!        write(6,*)r,F2(r)
!ENDDO

!CALL stdpara.readstd
!CALL spiral.init(500,12.d0,stdpara,2)
!CALL spiral.readw(2)
!CALL FindSpiral(spiral)
!DO i = 1, N
!        r = dr*dble(i)
!!       print *,spiral.r(i),&
!!       real(spiral.h1(i)),intplt(real(spiral.h1),spiral.r,r)
!        print *,r,intplt(real(spiral.h1),spiral.r,r)
!ENDDO

Y1(1) = dcmplx(0.d0,0.d0)
Y1(2) = dcmplx(1.d0,0.d0)
CALL rk4d1(0.d0,3.d0,20,F,F3,Y,Y1)
DO i = 1, 20
        print *,real(Y(1,i)),real(Y(2,i))
ENDDO

stop
CONTAINS

FUNCTION F(r)
IMPLICIT NONE
DOUBLE PRECISION                ::r
DOUBLE COMPLEX                  ::F
F = (0.d0,0.d0)
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

FUNCTION F3(r)
IMPLICIT NONE
DOUBLE PRECISION                ::r
DOUBLE COMPLEX                  ::F3
F3 = (1.d0,0.d0)*r
ENDFUNCTION

end      
