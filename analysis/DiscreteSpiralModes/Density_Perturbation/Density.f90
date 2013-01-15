PROGRAM spiral
USE PLOTTING
USE STELLARDISK,ToomreQ=>Q
IMPLICIT NONE
INTEGER                         ::i,j,k
CHARACTER(len=32)               ::arg
DOUBLE COMPLEX,ALLOCATABLE      ::u(:,:),h1(:),phi1r(:)
DOUBLE PRECISION                ::domain= 20.d0,dx,dy,r,th
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::potential(:,:)
DOUBLE PRECISION,ALLOCATABLE    ::force(:,:,:)
INTEGER,PARAMETER               ::n=500

!CALL getarg(1,arg)
!READ(arg,*)wr
!CALL getarg(2,arg)
!READ(arg,*)wi

 wr = 65.44d0
 wi = -0.71d0

ALLOCATE(u(3,4*n))
ALLOCATE(h1(4*n))
!find EigenFunction
CALL findu(u)
!find h1
CALL findh1(u,h1)
CALL refineh1(u,h1)
!!Find phi1 along r
ALLOCATE(phi1r(2*n))
CALL FindPhi1(phi1r)

dx = domain/dble(n)
dy = domain/dble(n)

ALLOCATE(density(2*n,2*n))
ALLOCATE(xcoord(2*n))
ALLOCATE(ycoord(2*n))

DO i = 1, n*2
        xcoord(i) =  dx*0.5d0 - domain + dble(i)*dx
        ycoord(i) =  dy*0.5d0 - domain + dble(i)*dx
ENDDO

DO i = 1, n*2
DO j = 1, n*2
        r = sqrt(xcoord(i)**2+ycoord(j)**2)
        th = atan2(ycoord(j),xcoord(i))
        density(i,j) = sigma1(r,th)
!       density(i,j) = sigma0(r)
ENDDO
ENDDO


open(10,file='r-dep.dat')
DO i = 2, n*4,2
         r = real(u(1,i))
        write(10,'(5(1XE15.6))')real(u(1,i)),real(u(2,i)),real(h1(i))/snsd(r)**2*sigma0(r),real(phi1r(i/2)),real(h1(i))
!        write(10,'(4(1XE15.6))')real(u(1,i)),real(phi1r(i/2))
enddo
close(10)

!!Find 2d Potential
ALLOCATE(potential(2*n,2*n))
DO i = 1, n*2
DO j = 1, n*2
        r = sqrt(xcoord(i)**2+ycoord(j)**2)
        th = atan2(ycoord(j),xcoord(i))
        potential(i,j) = phi1(r,th)
ENDDO
ENDDO

!!Find Force
ALLOCATE(force(n/4,n/4,2))
!DO i = 1,n/4
!DO j = 1, n/4
!        r = sqrt(xcoord(i*8)**2+ycoord(j*8)**2)
!        th = atan2(ycoord(j*8),xcoord(i*8))
!        call FindForce(force(i,j,:),r,th)
!ENDDO
!ENDDO

CALL plot2d(density,potential,force,n,domain)
DEALLOCATE(potential)
DEALLOCATE(xcoord)
DEALLOCATE(ycoord)
DEALLOCATE(u)
DEALLOCATE(density)
STOP
CONTAINS

SUBROUTINE FindForce(force,r,th)
IMPLICIT NONE
DOUBLE PRECISION                ::force(2),r,th
DOUBLE COMPLEX                  ::hh1,phir
DOUBLE PRECISION                ::runit(2),thunit(2)
INTEGER                         ::i
!interploting u at non-grid point r
do i = 1, n*2
        if(real(u(1,2*i)).gt.r)then
                exit
        endif
enddo
phir =(r-u(1,2*i-2))/(u(1,2*i)-u(1,2*i-2))*(phi1r(i)-phi1r(i-1)) + phi1r(i-1)
i = i*2
hh1 =(r-u(1,i-1))/(u(1,i)-u(1,i-1))*(h1(i)-h1(i-1)) + h1(i-1)

runit = (/cos(th),sin(th)/)
thunit = (/-sin(th),cos(th)/)

force = real(dsimplifiedPoisson(r,phir,hh1)*exp((0.d0,-2.d0)*th)*runit &
      + (0.d0,-2.d0)*phir/r*exp((0.d0,-2.d0)*th)*thunit)

ENDSUBROUTINE

FUNCTION p(r)
IMPLICIT NONE
DOUBLE COMPLEX                  ::p
DOUBLE PRECISION                ::r
p = (0.d0,0.d0)
ENDFUNCTION

FUNCTION q(r)
IMPLICIT NONE
DOUBLE COMPLEX                  ::q
DOUBLE PRECISION                ::r
q = k3sqrt(r)
!q = +1.d0
ENDFUNCTION

SUBROUTINE findu(u)
IMPLICIT NONE
DOUBLE COMPLEX                  ::u(:,:)
DOUBLE COMPLEX                  ::ui(3)
DOUBLE PRECISION                ::a,b


a = 0.00001d0
b = 2.d0*domain
ui = (/a,1.d0,0.d0/)
CALL rk4(a,b,4*N,p,q,p,u,ui)

ENDSUBROUTINE

SUBROUTINE refineh1(u,h1)
IMPLICIT NONE
DOUBLE COMPLEX                  ::u(:,:),h1(:)
INTEGER                         ::i,f
DO i =1, size(u,2)
        if(abs(u(1,i)).gt.9.5d0)then
                f = i
                exit
        endif
enddo

DO i = f+1, size(u,2)
        h1(i) = h1(i)*exp((-abs(u(1,i))+9.5d0)/0.5d0)
enddo

ENDSUBROUTINE

SUBROUTINE findh1(u,h1)
IMPLICIT NONE
DOUBLE COMPLEX                  ::u(:,:),h1(:)
DOUBLE COMPLEX                  ::rad
DOUBLE PRECISION                ::h,r
INTEGER                         ::i,j

do i = 1, 4*N 
!find h1 by the interploted u
        r   = u(1,i)
        rad = sqrt(kappa(r)**2*(1.d0-nu(r)**2)/sigma0(r)/r)
        h1(i) = u(2,i)*rad*exp(-0.5*(0.d0,1.d0)*ExpPart(r))
enddo

ENDSUBROUTINE

Function ExpPart(r)
IMPLICIT NONE
DOUBLE PRECISION                ::ExpPart,r
DOUBLE PRECISION                ::err = 10.d-10
INTEGER                         ::IERR
CALL DGAUS8(Sigma,0.d0,r,ERR,ExpPart,IERR)

ENDFUNCTION

function Sigma(r)
!This is NOT related to density
IMPLICIT NONE
DOUBLE PRECISION                ::Sigma,r
Sigma = 2.d0*pi*G*sigma0(r)/snsd(r)**2
ENDFUNCTION

FUNCTION sigma1(r,th)
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
DOUBLE PRECISION                ::sigma1
DOUBLE PRECISION,INTENT(IN)     ::r,th
DOUBLE COMPLEX                  ::uu,hh1
DOUBLE PRECISION                ::rad
INTEGER                         ::i,j,k,l
!interploting u at non-grid point r
do i = 1, n*3
        if(real(u(1,i)).gt.r)then
                exit
        endif
enddo
!uu = (r-u(1,i-1))/(u(1,i)-u(1,i-1))*(u(2,i)-u(2,i-1)) + u(2,i-1)
hh1 =(r-u(1,i-1))/(u(1,i)-u(1,i-1))*(h1(i)-h1(i-1)) + h1(i-1)


!find density
sigma1 = real(hh1*sigma0(r)/snsd(r)**2*exp(-2.d0*th*(0.d0,1.d0)))


ENDFUNCTION

FUNCTION phi1(r,th)
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
DOUBLE PRECISION                ::phi1
DOUBLE PRECISION,INTENT(IN)     ::r,th
DOUBLE COMPLEX                  ::uu,hh1
DOUBLE PRECISION                ::rad
INTEGER                         ::i,j,k,l
!interploting u at non-grid point r
do i = 1, n*2
        if(real(u(1,2*i)).gt.r)then
                exit
        endif
enddo
!uu = (r-u(1,i-1))/(u(1,i)-u(1,i-1))*(u(2,i)-u(2,i-1)) + u(2,i-1)
hh1 =(r-u(1,2*i-2))/(u(1,2*i)-u(1,2*i-2))*(phi1r(i)-phi1r(i-1)) + phi1r(i-1)


!find potential
phi1 = real(hh1*exp(-2.d0*th*(0.d0,1.d0)))


ENDFUNCTION

SUBROUTINE FindPhi1(phi1)
IMPLICIT NONE
DOUBLE COMPLEX                  ::phi1(:),k(4)
DOUBLE PRECISION                ::r,h
INTEGER                         ::i,j,l,nn

!!Solve ODE of phi from density by RK4
h = u(1,3)-u(1,1)
phi1(1) = (0.d0,0.d0)
DO i = 2,4*N-2,2
        r = u(1,i-1)
        k(1) = h*dsimplifiedPoisson(r           ,phi1(i/2)                ,h1(i-1))
        k(2) = h*dsimplifiedPoisson(r+h/2.d0    ,phi1(i/2)+k(1)/2.d0      ,h1(i))
        k(3) = h*dsimplifiedPoisson(r+h/2.d0    ,phi1(i/2)+k(2)/2.d0      ,h1(i))
        k(4) = h*dsimplifiedPoisson(r+h         ,phi1(i/2)+k(3)           ,h1(i+1))
        phi1(i/2+1) = phi1(i/2) + (k(1)+2.d0*k(2)+2.d0*k(3)+k(4))/6.d0
ENDDO


ENDSUBROUTINE

FUNCTION dsimplifiedPoisson(r,phi,h)
IMPLICIT NONE
DOUBLE COMPLEX                  ::dsimplifiedPoisson
DOUBLE COMPLEX,INTENT(IN)       ::phi,h
DOUBLE PRECISION                ::r

dsimplifiedPoisson = -phi/(2.d0*r)+(0.d0,1.d0)*cmplx(Sigma(r),0)*h
!dsimplifiedPoisson = dsimplifiedPoisson + 3.75d0/(0.d0,1.d0)/Sigma(r)/r**2*phi

!!test case 
!!dsimplifiedPoisson =  -phi*r + (0.d0,1.d0)*r
!!bnd condition is phi = 1+i at r=0
!!solution is
!!phi = exp(-r**2/2)+ i
!!dsimplifiedPoisson = -phi*r
ENDFUNCTION

END PROGRAM

