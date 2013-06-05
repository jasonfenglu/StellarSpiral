subroutine init2d(q_loc)
use common_params
use GALAXY,only:gasdensity,RC,spiral
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
integer::i,j
! implement your initial condition here
double precision:: r_loc,vx,vy,rho,CV
double precision:: pspd
pspd = real(spiral.w)/2.d0
q_loc = 0.d0
   do j=1,ncell_loc(2)
     do i=1,ncell_loc(1)
          r_loc = dsqrt(x_loc(i)**2+y_loc(j)**2)
          CV = RC(r_loc) - pspd*r_loc
          vx = -CV*y_loc(j)/r_loc
          vy =  CV*x_loc(i)/r_loc
          rho = gasdensity(r_loc,0.d0)
          q_loc(i,j,1) = rho
          q_loc(i,j,2) = vx*rho
          q_loc(i,j,3) = vy*rho
     enddo
   enddo

end subroutine