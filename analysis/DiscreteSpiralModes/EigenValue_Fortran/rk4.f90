!   BEGIN PPROLOGUE RK4
!   PURPOSE   SOVLE SECONND ORDER ODE BY RK45
!   LIBRARY   
!   CATEGORY  ODE
!   TYPE      RK45
!   KEYWORDS  ODE, RK45
!   AUTHOR    Chien-Chang Feng
!   DESCRIPTION
!
!   Solve ODE like:
!       u''[r] + p(r) u' + q(r) u = f(r)
!
!   Implicit Runge-Kutta methods
!       y_(n+1) = y_n +h* sum(b_i kj)
!       kj =  f(t_n + c_i*h, y_n + h* sum(a_ij kj))
!
!   Complex Ready
!
!   REFERENCES  RKF_ABM.PDF(ON WEB)
!   ROUTINES CALLED  
!   REVISION HISTORY  (YYMMDD)
!   121122   DATE WRITTEN
!***END PROLOGUE  RK4

MODULE RK
IMPLICIT NONE
CONTAINS
SUBROUTINE rk45(ri,rf,N,p,q,s,u,ui)
IMPLICIT NONE
DOUBLE COMPLEX          ::u(3,N),ui(3),uu(3)
DOUBLE COMPLEX          ::K(3,7)
DOUBLE COMPLEX,EXTERNAL ::p,q,s
DOUBLE PRECISION        ::a(6,7)
DOUBLE PRECISION        ::C(6)
DOUBLE PRECISION        ::b(7)
DOUBLE PRECISION        ::r
DOUBLE PRECISION        ::ri,rf,h
INTEGER                 ::N
INTEGER                 ::I,J,L


!init boundary
h = (rf-ri)/REAL(N)
r = ri

!init RKF Tableau
a = 0.d0

a(:,1) = C

a(2,2) = 0.25d0

a(3,2) = 3.d0/32.d0
a(3,3) = 9.d0/32.d0

a(4,2) = 1932.d0/2197.d0
a(4,3) = -7200.d0/2197.d0
a(4,4) = 7296.d0/2197.d0

a(5,2) = 439.d0/216.d0
a(5,3) = -8.d0
a(5,4) = 3680.d0/513.d0
a(5,5) = -845.d0/4104.d0

a(6,2) = -8.d0/27.d9
a(6,3) = 2.d0
a(6,4) = -3544.d0/2565.d0
a(6,5) = 1859.d0/4104.d0
a(6,6) = -11.d0/40.d0

b = (/0.d0,25.d0/216.d0,0.d0,1408.d0/2565.d0,2197.d0/4104.d0,-0.2d0,0.d0/)
C = (/0.d0,0.25d0,3.d0/8.d0,12.d0/13.d0,1.d0,0.5d0/)

u(:,1) = ui
!iteration u
do l = 2,n
 !iterating k
  k = 0.d0
  do i = 1,7 
          uu = u(:,l)
          do j = 1,7
                  uu = uu + a(i,j)*k(:,j)
          enddo
          k(:,i) = f(u(:,l-1) + h*uu,p,q,s)
  enddo

  uu = 0.d0
  do j = 1,7
          uu = uu + b(j)*k(:,j)
  enddo

  u(:,l) = u(:,l-1) + uu*h

enddo


write(*,*)REAL(u)


contains
FUNCTION f(ui,p,q,s)
IMPLICIT NONE
DOUBLE COMPLEX          ::f(3)
DOUBLE COMPLEX          ::ui(3)
DOUBLE COMPLEX,EXTERNAL::p,q,s
f(1) = 1.d0
f(2) = ui(3)
f(3) = s(r)-p(r)*ui(3)-q(r)*ui(2)
ENDFUNCTION

ENDSUBROUTINE

SUBROUTINE rk4(ri,rf,N,p,q,s,u,ui)
IMPLICIT NONE
include 'omp_lib.h'
DOUBLE COMPLEX          ::u(3,N),ui(3),uu(3)
DOUBLE COMPLEX          ::K(3,4)
DOUBLE COMPLEX,EXTERNAL ::p,q,s
DOUBLE PRECISION        ::a(4,5)
DOUBLE PRECISION        ::C(4)
DOUBLE PRECISION        ::b(4)
DOUBLE PRECISION        ::r
DOUBLE PRECISION        ::ri,rf,h
INTEGER                 ::N
INTEGER                 ::I,J,L

!init boundary
h = (rf-ri)/REAL(N)
r = ri

!init RKF Tableau

a = 0.d0
a(2,1) = 0.5d0
a(2,2) = 0.5d0
a(3,1) = 0.5d0
a(3,3) = 0.5d0
a(4,1) = 1.d0
a(4,4) = 1.d0

b = (/1.d0/6.d0,1.d0/3.d0,1.d0/3.d0,1.d0/6.d0/)

C = (/0.d0,0.5d0,0.5d0,1.d0/)


u(:,1) = ui

!iteration u
do l = 2,n
 !iterating k
 k = 0.d0
 !find k_i
 do i = 1,4
         uu = u(:,l)
         do j =  1,5
                uu = uu + a(i,j)*k(:,j)
         enddo
         k(:,i) = f(u(:,l-1) + h*uu,p,q,s)
 enddo

 uu = 0.d0
 do j = 1,4
        uu = uu + b(j)*k(:,j)
 enddo

 u(:,l) = u(:,l-1) + h*uu
enddo




contains
RECURSIVE FUNCTION f(ui,p,q,s) result(ans)
IMPLICIT NONE
DOUBLE COMPLEX          ::ans(3)
DOUBLE COMPLEX          ::ui(3)
DOUBLE COMPLEX,EXTERNAL::p,q,s
DOUBLE PRECISION        ::r
r = ui(1)
ans(1) = 1.d0
ans(2) = ui(3)
ans(3) = s(r)-p(r)*ui(3)-q(r)*ui(2)
ENDFUNCTION

ENDSUBROUTINE



!do i = 1,n-1
!        ui     = u(:,i)
!        ui(1)  = r
!        u(1,i) = r
!
!        K(:,1) = h*f(ui)
!
!        w = 0.d0
!        w(1) = 0.25d0*h
!        w(2) = 0.25d0
!        K(:,2) = h*f(ui +sum( K*DCMPLX(spread(w,1,3)),2))
!
!        w = 0.d0
!        w(1) = 3.d0/8.d0*h
!        w(2) = 3.d0/32.d0
!        w(3) = 9.d0/32.d0
!        K(:,3) = h*f(ui +sum( K*DCMPLX(spread(w,1,3)),2))
!
!        w = 0.d0
!        w(1) = 12.d0/13.d0*h
!        w(2) = 1932.d0/2197.d0
!        w(3) = -7200.d0/2197.d0
!        w(4) = 7296.d0/2197.d0
!        K(:,4) = h*f(ui +sum( K*DCMPLX(spread(w,1,3)),2))
!
!        w = 0.d0
!        w(1) = 1.d0*h
!        w(2) = 439.d0/216.d0
!        w(3) = -8.d0
!        w(4) = 3680.d0/513.d0
!        w(5) = -845.d0/4104.d0
!        K(:,5) = h*f(ui +sum( K*DCMPLX(spread(w,1,3)),2))
!
!        w = 0.d0
!        w(1) = 0.5d0*h
!        w(2) = -8.d0/27.d0
!        w(3) = 2.d0
!        w(4) = -3544.d0/2565.d0
!        w(5) = 1859.d0/4104.d0
!        w(6) = -11.d0/40.d0
!        K(:,6) = h*f(ui +sum( K*DCMPLX(spread(w,1,3)),2))
!
!        w = 0.d0
!        w(2) = 25.d0/216.d0
!        w(4) = 1408.d0/2565.d0
!        w(5) = 2197.d0/4104.d0
!        w(6) = -1.d0/5.d0
!        w = w *h
!        u(:,i+1) = u(:,i) +sum( K*DCMPLX(spread(w,1,3)),2)
!
!        r = r + h
!
!enddo

ENDMODULE RK
