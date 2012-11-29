PROGRAM find_one
IMPLICIT NONE
INTEGER                         ::i
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::wr,wi
DOUBLE PRECISION                ::error

CALL getarg(1,arg)
READ(arg,*)wr
CALL getarg(2,arg)
READ(arg,*)wi

write(*,*)wr,wi
do i = 1, 1000
        CALL update(wr,wi,error)
        write(*,*)wr,wi,error
enddo

STOP

CONTAINS

FUNCTION p(r)
IMPLICIT NONE
DOUBLE COMPLEX                  ::p
DOUBLE PRECISION                ::r
p = (0.d0,0.d0)
ENDFUNCTION

FUNCTION q(r)
IMPLICIT NONE
DOUBLE COMPLEX                  ::q
DOUBLE COMPLEX,EXTERNAL         ::k3sqrt
DOUBLE PRECISION                ::r
q = k3sqrt(r,wr,wi)
ENDFUNCTION

SUBROUTINE update(wr,wi,err)
IMPLICIT NONE
DOUBLE COMPLEX,ALLOCATABLE      ::u(:,:),ui(:)
DOUBLE COMPLEX,EXTERNAL         ::error,k3sqrt
DOUBLE PRECISION,EXTERNAL       ::find_b,KappaOverASqr,mo,kappa
DOUBLE PRECISION                ::a,b
DOUBLE PRECISION                ::h
DOUBLE PRECISION                ::wr,wi
DOUBLE PRECISION,optional       ::err
DOUBLE COMPLEX                  ::nu,nusqr
INTEGER                         ::N=100
INTEGER                         ::i
DOUBLE PRECISION                ::diff


ALLOCATE(u(3,N))
ALLOCATE(ui(3))
a = 0.d0
b = find_b(wr)
ui = (/a,1.d0,0.d0/)
CALL rk4(a,b,N,p,q,p,u,ui)
err = ABS(u(3,N)/u(2,N)-error(b,wr,wi))

h = 10.d-15
nusqr= u(3,N)/u(2,N)+ &
       0.5d0/sqrt(k3sqrt(b,wr,wi)) &
       *(sqrt(k3sqrt(b+h,wr,wi))-sqrt(k3sqrt(b-h,wr,wi)))/(2.d0*h)

nusqr = -nusqr**2
nu    = nusqr/KappaOverASqr(b) + 1.d0 -( 1.d0 + 0.75d0*dexp(-b**2))**-2
nu    = sqrt(nusqr)*kappa(b)+mo(b)

diff = real(nu) - wr
diff = abs(diff)/diff*sqrt(diff)
wr = wr+(real(nu) - wr)*10.d-5

diff = aimag(nu) - wi
diff = abs(diff)/diff*sqrt(diff)
wi = wi+(aimag(nu) - wi)*10.d-5

DEALLOCATE(u)
DEALLOCATE(ui)
ENDSUBROUTINE


END PROGRAM

