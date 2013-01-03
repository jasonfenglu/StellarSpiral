        module plotting
        CONTAINS
        SUBROUTINE plot2d(F,F2,n,domain)
        IMPLICIT NONE
        DOUBLE PRECISION,ALLOCATABLE,INTENT(IN) ::F(:,:),F2(:,:)!plotting data
        DOUBLE PRECISION                        ::domain!plot range
        REAL                                    ::TR(6) !plot geometry
        REAL                                    ::vmax,vmin
        REAL                                    ::BRIGHT,CONTRA
        INTEGER                                 ::m,n   !dimentsion
        INTEGER                                 ::PGBEG
        REAL                                    ::dx,dy


        IF (PGBEG(0,'/png',1,1) .NE. 1) STOP
        CALL PGSVP(0.0,0.95,0.0,0.95)

        m = n
        dx = real(domain)/real(n)
        dy = real(domain)/real(m)

        TR(3) = 0.
        TR(5) = 0.
        TR(2) = 1.d0/(2.d0*dble(n)-1.)*(2.d0*domain-dx)
        TR(1) = TR(2)*dble(n)-2.*domain-dx
        TR(6) = 1.d0/(2.d0*dble(m)-1.)*(2.d0*domain-dy)
        TR(4) = TR(6)*dble(m)-2.*domain-dy


!       TR = (/0.,1.,0.,0.,0.,1./)
        vmax = real(MAXVAL(F(:,:)))
        vmin = real(MINVAL(F(:,:)))
!       write(*,*)vmax,vmin

!       vmax = 4.d6
!       vmin = -4.d6

        BRIGHT = 0.5
        CONTRA = -0.9
        CALL PALETT(2,CONTRA,Bright)

!        CALL PGBBUF
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        CALL PGIMAG(REAL(F(:,:)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmax, vmin, 'pixel value')
        CALL PGSCH(1.0)
!       CALL PGEBUF


!       vmax = real(MAXVAL(F2(:,:)))
!       vmin = real(MINVAL(F2(:,:)))
!       write(*,*)vmax,vmin

!       vmax = 10.d4
!       vmin = -10.d4
!       CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
!       CALL PGIMAG(REAL(F2),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
!       CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, 'pixel value')
!       CALL PGSCH(1.0)

        CALL PGCLOS

        ENDSUBROUTINE


      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
!-----------------------------------------------------------------------
! Set a "palette" of colors in the range of color indices used by
! PGIMAG.
!-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
!
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
!
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
!
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
!
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
               0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
               0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
               0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
      IF (TYPE.EQ.1) THEN
!        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
!        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
!        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
!        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
!        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      ENDSUBROUTINE

        endmodule

PROGRAM spiral
USE PLOTTING
USE STELLARDISK,ToomreQ=>Q
IMPLICIT NONE
INTEGER                         ::i,j,k
CHARACTER(len=32)               ::arg
DOUBLE COMPLEX,ALLOCATABLE      ::u(:,:),h1(:),phi1r(:)
DOUBLE PRECISION                ::domain= 12.d0,dx,dy,r,th
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::potential(:,:)
INTEGER,PARAMETER               ::n=500

!CALL getarg(1,arg)
!READ(arg,*)wr
!CALL getarg(2,arg)
!READ(arg,*)wi

 wr = 59.218d0
 wi = -0.855d0
!wr = 47.393d0
!wi = -0.533d0
!wr = 39.500d0
!wi = -0.400d0

ALLOCATE(u(3,4*n))
ALLOCATE(h1(4*n))
!find EigenFunction
CALL findu(u)
!find h1
CALL findh1(u,h1)
!find phi
!DO i = 1, n*3
! write(*,*)real(u(1,i)),real(h1(i))
!enddo
!stop

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

!CALL TEST(DENSITY)

CALL plot2d(density,density,n,domain)
ALLOCATE(phi1r(2*n))
CALL FindPhi1(phi1r)

open(10,file='r-dep.dat')
DO i = 2, n*4,2
         write(10,'(4(1XE15.6))')real(u(1,i)),real(u(2,i)),real(h1(i)),real(phi1r(i/2))
enddo
close(10)


ALLOCATE(potential(2*n,2*n))
DO i = 1, n*2
DO j = 1, n*2
        r = sqrt(xcoord(i)**2+ycoord(j)**2)
        th = atan2(ycoord(j),xcoord(i))
        potential(i,j) = phi1(r,th)
ENDDO
ENDDO
 CALL plot2d(potential,potential,n,domain)


!DEALLOCATE(potential)
DEALLOCATE(xcoord)
DEALLOCATE(ycoord)

DEALLOCATE(u)
DEALLOCATE(density)
STOP
CONTAINS

SUBROUTINE FindPotential(density,potential)
use,intrinsic                   ::iso_c_binding
IMPLICIT NONE
include                         'fftw3.f03'
DOUBLE PRECISION                ::density(:,:),potential(:,:)

type(C_PTR)                     ::plan
complex(C_DOUBLE_COMPLEX),ALLOCATABLE::in(:,:),out(:,:)
INTEGER                         ::n
INTEGER                         ::i,j

n = size(density,1)*2
ALLOCATE(in(n,n))
ALLOCATE(out(n,n))

in = 0.d0
do i = 1, n/2
do j = 1, n/2
        in(i,j) = DCMPLX(density(i,j))
enddo
enddo

plan = fftw_plan_dft_2d(n,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
CALL fftw_execute_dft(plan,in,out)
call fftw_destroy_plan(plan)


 in = (0.d0,0.d0)
 do i = 1, n
 do j = 1, n
         in(i,j) = out(i,j)/(-pi)/(dcmplx(i)**2+dcmplx(j)**2)*4.d0*g
 enddo
 enddo
 
 plan = fftw_plan_dft_2d(n,n,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
 CALL fftw_execute_dft(plan,in,out)
 call fftw_destroy_plan(plan)

do i = 1, n/2
do j = 1, n/2
        potential(i,j) = -1.d0*dble(out(i,j))/dble(n**2)
enddo
enddo

DEALLOCATE(in)
DEALLOCATE(out)
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

SUBROUTINE findh1(u,h1)
IMPLICIT NONE
DOUBLE COMPLEX                  ::u(:,:),h1(:)
DOUBLE PRECISION                ::rad,h,r
INTEGER                         ::i,j

do i = 1, 4*N 
!find h1 by the interploted u
        r   = u(1,i)
        rad = sqrt(kappa(r)**2*(1-nu(r)**2)/sigma0(r)/r)
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


!find density
phi1 = real(hh1*exp(-2.d0*th*(0.d0,1.d0)))


ENDFUNCTION

SUBROUTINE FindPhi1(phi1)
IMPLICIT NONE
DOUBLE COMPLEX                  ::phi1(:),k(4)
DOUBLE PRECISION                ::r,h
INTEGER                         ::i,j,l,nn

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
dsimplifiedPoisson = -phi/(2.d0*r)+(0.d0,1.d0)*cmplx(Sigma(r))*h

!!test case 
!!dsimplifiedPoisson =  -phi*r + (0.d0,1.d0)*r
!!bnd condition is phi = 1+i at r=0
!!solution is
!!phi = exp(r**2/2)+ i
ENDFUNCTION

SUBROUTINE test(density)
IMPLICIT NONE
DOUBLE PRECISION                ::density(:,:)
DOUBLE PRECISION                ::x,y,a,b
INTEGER                        ::t,i,j

a = 1.d0
b = 1.d0
density=0.d0

do t = 1,1000
        x = a*dexp(b*dble(t)/10.d1)*dcos(dble(t/10.d1))
        y = a*dexp(b*dble(t)/10.d1)*dsin(dble(t/10.d1))
        if(abs(x).lt.domain .and. abs(y).lt.domain)then
                do i = 1,2*n
                       if(xcoord(i).gt.x)then
                                exit
                       endif
                enddo
                do j = 1,2*n
                       if(ycoord(j).gt.y)then
                                exit
                       endif
                enddo

                density(i,j) = 10.d0
        write(*,*)x,y
        endif
enddo

!!!test plot
!DO i = 1,2*n
!DO j = 1,2*n
!        
!        x = xcoord(i)
!        y = ycoord(j)
!        density(i,j) = 1.d0*(x-3.d0)*y
!
!ENDDO
!ENDDO
!!!

ENDSUBROUTINE



END PROGRAM

