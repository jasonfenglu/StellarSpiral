subroutine force2d(q_loc,fx,fy)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::fx(1:ncell_loc(1),1:ncell_loc(2)),fy(1:ncell_loc(1),1:ncell_loc(2))
#ifdef GRAVITY
double precision,dimension(:,:),allocatable::den
#endif
double precision::r_loc,CV,CV2,Orsq
double precision::pspd,bfam,a1,p0,rsq,p1,psi,dpsi,picle
double precision::cs2s,sn2s,cs2t,sn2t,cs2st,sn2st
double precision::px,py
double precision::sndsq,density_distr,rho,vx,vy
integer::i,j
double precision::A,B,maxrho
pspd=21.6385d0
maxrho=15.d0
A=-0.249125d0
B=0.552313d0


! implement your external force here !!
  do j=1, ncell_loc(2)
    do i=1, ncell_loc(1)
       r_loc = dsqrt(x_loc(i)**2.d0+y_loc(j)**2.d0)
       rho = q_loc(i,j,1)
       vx = q_loc(i,j,2)/rho
       vy = q_loc(i,j,3)/rho
       !CV = 224.d0*dsqrt(r_loc/(r_loc+0.01))
       CV = 549.959d0*(r_loc/(r_loc**B+r_loc**(1-A)))
       sndsq= snd**2.d0
!       Orsq= (CV/r_loc)**2.d0 + 2.d0*sndsq*dexp(-(r_loc/6.98986d0)**2)/(6.98986d0**2.d0)
       Orsq= (CV/r_loc)**2.d0 + 2.d0*sndsq/(14.311d0**2.d0)
       !Orsq= (CV/r_loc)**2.d0
       fx(i,j)=-x_loc(i)*Orsq + (pspd**2.d0)*x_loc(i) + 2.d0*pspd*vy
       fy(i,j)=-y_loc(j)*Orsq + (pspd**2.d0)*y_loc(j) - 2.d0*pspd*vx
    enddo
  enddo

#ifdef GRAVITY
allocate(den(1:ncell_loc(1),1:ncell_loc(2)))
  do j=1, ncell_loc(2)
    do i=1, ncell_loc(1)
      r_loc= dsqrt(x_loc(i)**2.d0+y_loc(j)**2.d0)
      density_distr= maxrho*dexp(-(r_loc/14.311d0)**2.d0)
      !density_distr= 1.53d0
      den(i,j)= q_loc(i,j,1)-density_distr
     ! den= q_loc(1:ncell_loc(1),1:ncell_loc(2),1)
    enddo
  enddo
  call selfgravity2d(den)
  fx = fx+sgx
  fy = fy+sgy 
deallocate(den)
#endif

#ifdef BAR
!Add bar potential(type 2 bar)
 !pspd = 25.d0
 a1 = 2.2d0
 bfam = -255297.d0 !21% a1=2.2
 picle = 3.1415926535897932d0
! p0 = bfam*dmin1(t/tend*5.d0,1.0d0)
! p0 = bfam*dmin1(0.5d0*(1.0d0-dcos(picle*t/tend*5.d0)),1.0d0)
if (t .le. tend/5.d0) then
 p0 = bfam*0.5d0*(1.0d0-dcos(picle*t/tend*5.d0))
else
 p0 = bfam
endif
 

do j=1, ncell_loc(2)
   do i=1, ncell_loc(1)
      rsq = x_loc(i)**2.d0+y_loc(j)**2.d0
      r_loc = dsqrt(rsq)
       p1   = rsq/(a1**2+rsq)**2
      psi   = p0*p1
      dpsi  = p0*2.d0*r_loc*(a1**2-rsq)/(a1**2+rsq)**3
      cs2s  = (x_loc(i)**2-y_loc(j)**2)/rsq
      sn2s  = 2.d0*y_loc(j)*x_loc(i)/rsq
      !cs2t  = cos(2.d0*pspd*t)
      !sn2t  = sin(2.d0*pspd*t)
      !cs2st = cs2s*cs2t+sn2s*sn2t
      !sn2st = sn2s*cs2t-cs2s*sn2t
      !px    = -((dpsi*cs2st)*x_loc(i)/r_loc+2.d0*psi*sn2st*y_loc(j)/rsq)
      !py    = -((dpsi*cs2st)*y_loc(j)/r_loc-2.d0*psi*sn2st*x_loc(i)/rsq)
      px    = -((dpsi*cs2s)*x_loc(i)/r_loc+2.d0*psi*sn2s*y_loc(j)/rsq)
      py    = -((dpsi*cs2s)*y_loc(j)/r_loc-2.d0*psi*sn2s*x_loc(i)/rsq)
    fx(i,j) = fx(i,j) + px
    fy(i,j) = fy(i,j) + py
   enddo
 enddo 
#endif

end subroutine

