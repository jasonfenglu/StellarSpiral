subroutine force2d(q_loc,fx,fy)
use common_params
use simcontroll
use GALAXY, only:RC,stellarforce,spiral
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::fx(1:ncell_loc(1),1:ncell_loc(2)),fy(1:ncell_loc(1),1:ncell_loc(2))
#ifdef GRAVITY
double precision,dimension(:,:),allocatable::den
#endif
integer::i,j
!!gas
double precision::r_loc,CV,Orsq
!!bar
!double precision::pspd,bfam,a1,p0,rsq,p1,psi,dpsi
!double precision::cs2s,sn2s,cs2t,sn2t,cs2st,sn2st
!double precision::px,py
!!rotating frame
DOUBLE PRECISION::pspd
DOUBLE PRECISION::rho,vx,vy

pspd = real(spiral.w)/2.d0
! implement your external force here !!
! gas pressure
  do j=1, ncell_loc(2)
    do i=1, ncell_loc(1)
       rho = q_loc(i,j,1)
       vx = q_loc(i,j,2)/rho
       vy = q_loc(i,j,3)/rho
       r_loc = dsqrt(x_loc(i)**2.d0+y_loc(j)**2.d0)
       CV = RC(r_loc)
       Orsq= (CV/r_loc)**2.d0
       !fx(i,j)=-x_loc(i)*Orsq
       !fy(i,j)=-y_loc(j)*Orsq
       fx(i,j)=-x_loc(i)*Orsq + (pspd**2.d0)*x_loc(i) + 2.d0*pspd*vy
       fy(i,j)=-y_loc(j)*Orsq + (pspd**2.d0)*y_loc(j) - 2.d0*pspd*vx
    enddo
  enddo


!!old spiral potential, not using 
!do j=1, ncell_loc(2)
!  do i=1, ncell_loc(1)
!     r_loc = dsqrt(x_loc(i)**2.d0+y_loc(j)**2.d0)
!     th    = atan2(y_loc(j),x_loc(i))
!      CALL FindForce(fspi,r_loc,th-wr*t/2.d0)
!      !force is 1%
!      fspi = fspi * dmin1(t/tend*dble(simcon%ncir),1.d0) *2.5d0
!      fx(i,j)= fspi(1) + fx(i,j)
!      fy(i,j)= fspi(2) + fy(i,j)
!      force(i,j,:) = fspi(:)
!  enddo
!enddo

#ifdef STELLARSPIRAL
if(myid.eq.0)then
        write(6,*)achar(27)//'[95m applying stellar force'//achar(27)//'[0m'
ENDIF
fx = fx + stellarforce.sgx * dmin1(t/tend*dble(simcon.ncir),1.d0)
fy = fy + stellarforce.sgy * dmin1(t/tend*dble(simcon.ncir),1.d0)
if(myid.eq.0)then
        write(6,*)achar(27)//'[95m stellar complete'//achar(27)//'[0m'
ENDIF
#endif

#ifdef GRAVITY
 allocate(den(1:ncell_loc(1),1:ncell_loc(2)))
 den= q_loc(1:ncell_loc(1),1:ncell_loc(2),1)-1.53d2
 !den= q_loc(1:ncell_loc(1),1:ncell_loc(2),1)
 call selfgravity2d(den)
 fx = fx+sgx
 fy = fy+sgy 
 deallocate(den)
#endif

!Add bar potential 
! pspd = 232.554d0
! a1 = 0.09d0
! bfam = -27.1301d0
! p0 = bfam*dmin1(t/tend*10.d0,1.d0)
! do j=1, ncell_loc(2)
!   do i=1, ncell_loc(1)
!      rsq = x_loc(i)**2.d0+y_loc(j)**2.d0
!      r_loc = dsqrt(rsq)
!       p1   = rsq/(a1**2+rsq)**2
!      psi   = p0*p1
!      dpsi  = p0*2.d0*r_loc*(a1**2-rsq)/(a1**2+rsq)**3
!      cs2s  = (x_loc(i)**2-y_loc(j)**2)/rsq
!      sn2s  = 2.d0*y_loc(j)*x_loc(i)/rsq
!      cs2t  = cos(2.d0*pspd*t)
!      sn2t  = sin(2.d0*pspd*t)
!      cs2st = cs2s*cs2t+sn2s*sn2t
!      sn2st = sn2s*cs2t-cs2s*sn2t
!      px    = -((dpsi*cs2st)*x_loc(i)/r_loc+2.d0*psi*sn2st*y_loc(j)/rsq)
!      py    = -((dpsi*cs2st)*y_loc(j)/r_loc-2.d0*psi*sn2st*x_loc(i)/rsq)
!    fx(i,j) = fx(i,j) + px
!    fy(i,j) = fy(i,j) + py
!   enddo
! enddo 

end subroutine
