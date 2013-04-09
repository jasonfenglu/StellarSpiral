program test
USE STELLARDISK
DOUBLE PRECISION        ::ri,rf,dr,r
INTEGER                 ::n,i

ri = 0.d0
rf = 1.d0
!rf = 10.d0
N  = 1000

dr = (rf-ri)/dble(n)

CALL INIT_STELLARDISK(n,rf)
DO i = 0, N
        r = dr*dble(i)
        write(6,*)r,real(curf(r))
ENDDO


stop
end      
