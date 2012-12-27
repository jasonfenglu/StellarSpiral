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
DOUBLE COMPLEX,ALLOCATABLE      ::u(:,:)
DOUBLE PRECISION                ::domain= 12.d0,dx,dy,r,th
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::potential(:,:)
INTEGER,PARAMETER               ::n=100

!CALL getarg(1,arg)
!READ(arg,*)wr
!CALL getarg(2,arg)
!READ(arg,*)wi

!wr = 59.218d0
!wi = -0.855d0
!wr = 47.393d0
!wi = -0.533d0
wr = 39.500d0
wi = -0.400d0

ALLOCATE(u(3,3*n))
!find EigenFunction
CALL findu(u)
DO i = 1, n*3
 write(*,*)real(u(1,i)),real(real(u(2,i)))
enddo

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
!$OMP PARALLEL SHARED(i,xcoord,ycoord,density) PRIVATE(j,r,th) DEFAULT(NONE)
!$OMP DO SCHEDULE(DYNAMIC)
DO j = 1, n*2
        r = sqrt(xcoord(i)**2+ycoord(j)**2)
        th = atan2(ycoord(j),xcoord(i))
        density(i,j) = sigma1(r,th)
!       density(i,j) = sigma0(r)
ENDDO
!$OMP END DO
!$OMP END PARALLEL
ENDDO

!CALL TEST(DENSITY)



ALLOCATE(potential(2*n,2*n))
!CALL FindPotential(density,potential)
!CALL plot2d(density,n,domain)
CALL plot2d(density,potential,n,domain)


!DEALLOCATE(potential)
DEALLOCATE(xcoord)
DEALLOCATE(ycoord)

DEALLOCATE(u)
DEALLOCATE(density)
STOP
CONTAINS

SUBROUTINE FindPotential(density,potential)
DOUBLE PRECISION                ::density(:,:),potential(:,:)
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
use rk
IMPLICIT NONE
DOUBLE COMPLEX                  ::u(:,:)
DOUBLE COMPLEX                  ::ui(3)
DOUBLE PRECISION                ::a,b


a = 0.0000001d0
b = 2.d0*domain
ui = (/a,1.d0,0.d0/)
CALL rk4(a,b,3*N,p,q,p,u,ui)

ENDSUBROUTINE

Function ExpPart(r)
IMPLICIT NONE
DOUBLE PRECISION                ::ExpPart,r
DOUBLE PRECISION                ::err = 10.d-10
INTEGER                         ::IERR
CALL DGAUS8(Sigma,0.d0,r,ERR,ExpPart,IERR)

ENDFUNCTION

function Sigma(r)
IMPLICIT NONE
DOUBLE PRECISION                ::Sigma,r
Sigma = 2.d0*pi*G*sigma0(r)/snsd(r)**2
ENDFUNCTION

FUNCTION sigma1(r,th)
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
DOUBLE PRECISION                ::sigma1
DOUBLE PRECISION,INTENT(IN)     ::r,th
DOUBLE COMPLEX                  ::uu,h1
DOUBLE PRECISION                ::rad
INTEGER                         ::i,j,k,l
!interploting u at non-grid point r
do i = 1, n*3
        if(real(u(1,i)).gt.r)then
                exit
        endif
enddo
uu = (r-u(1,i-1))/(u(1,i)-u(1,i-1))*(u(2,i)-u(2,i-1)) + u(2,i-1)

rad = sqrt(kappa(r)**2*(1-nu(r)**2)/sigma0(r)/r)

h1 = uu*rad*exp(-0.5*(0.d0,1.d0)*ExpPart(r)-2.d0*th*(0.d0,1.d0))

sigma1 = real(h1)*sigma0(r)/snsd(r)

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

