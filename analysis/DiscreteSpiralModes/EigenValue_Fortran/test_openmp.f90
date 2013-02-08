program test
USE STELLARDISK
USE OMP_LIB
IMPLICIT NONE
DOUBLE PRECISION                ::r,ans(100)
INTEGER                         ::i

!$OMP PARALLEL 
CALL INIT_STELLARDISK(100,20.d0)
!$OMP BARRIER
!$OMP DO
do i = 1, 100
        wr = 63.d0
        wi = -0.5d0
        CALL FindSpiral
        print *,i
        ans(i) = error()
!       r = 1.d0/1000.d0*dble(i)+0.0001d0
!       print *,r,sigma(r)
enddo
!$OMP END DO
!$OMP END PARALLEL
do i = 1, 100
        print *,ans(i)
enddo

CALL ENDSTELLARDISK
endprogram
