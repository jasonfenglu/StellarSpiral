subroutine dt_loc3d(q_loc,dt_loc)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR)
double precision::dt_loc
double precision::den,vx,vy,vz,ene,pressure,wspd_max,wspd,vtot
#ifdef MHD
double precision::bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc,vsq,bsq,cfast,bmin
#endif
integer::i,j,k
! implement your local dt here

wspd_max=0.d0 !maxima wave speed
#ifdef MHD
bmin=1.d10
#endif
do k=1,ncell_loc(3)
  do j=1,ncell_loc(2)
    do i=1,ncell_loc(1)
      den=q_loc(i,j,k,1)
       vx=q_loc(i,j,k,2)/den
       vy=q_loc(i,j,k,3)/den
       vz=q_loc(i,j,k,4)/den
      ene=q_loc(i,j,k,5)
     vtot=dsqrt(vx**2+vy**2+vz**2)
#ifdef ISOTHERMAL
     wspd=vtot+snd
#endif
#ifdef ADIABATIC
#ifdef MHD
     vz=q_loc(i,j,k,4)/den
     bxl=q_loc(i,j,k,5)
     byl=q_loc(i,j,k,6)
     bzl=q_loc(i,j,k,7)
     ene=q_loc(i,j,k,8) 
     bxr=q_loc(i,j,k,9)
     byr=q_loc(i,j,k,10)
     bzr=q_loc(i,j,k,11)
     bxc=0.5d0*(bxl+bxr)
     byc=0.5d0*(byl+byr)
     bzc=0.5d0*(bzl+bzr) 
     vsq=vx**2.d0+vy**2.d0+vz**2.d0
     bsq=bxc**2.d0+byc**2.d0+bzc**2.d0
    pressure=(gam-1.d0)*(ene-0.5d0*den*vsq-0.5d0*bsq)
    bmin=dmin1(dmin1(dabs(bxc),dabs(byc)),dabs(bzc))   
    cfast=dsqrt((gam*pressure+bsq+dsqrt((gam*pressure+bsq)**2.d0-4.d0*gam*pressure*bmin**2.d0))/(2.d0*den))
     wspd=vtot+cfast
#else
     pressure=(gam-1.d0)*(ene-0.5d0*den*(vx**2+vy**2+bz**2))
     wspd=vtot+dsqrt(gam*pressure/den)
#endif
#endif
     wspd_max=dmax1(wspd_max,wspd)
    enddo
  enddo
enddo
dt_loc=CFL*1.d0/wspd_max*dx
end subroutine
