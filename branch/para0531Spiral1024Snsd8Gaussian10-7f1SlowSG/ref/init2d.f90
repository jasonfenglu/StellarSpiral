subroutine init2d(q_loc)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
integer::i,j
! implement your initial condition here
double precision:: r_loc,vx,vy,rho,CV,pspd,A,B,maxrho
q_loc = 0.d0
pspd = 21.6385d0
maxrho=15.d0
A=-0.249125d0
B=0.552313d0

   do j=1,ncell_loc(2)
     do i=1,ncell_loc(1)
          r_loc = dsqrt(x_loc(i)**2+y_loc(j)**2)
          !CV = 224.d0*dsqrt(r_loc/(r_loc+0.01))-r_loc*pspd
          CV = 549.959d0*r_loc/(r_loc**B+r_loc**(1-A))-r_loc*pspd
          vx = -CV*y_loc(j)/r_loc
          vy =  CV*x_loc(i)/r_loc
          rho = maxrho*dexp(-(r_loc/14.311d0)**2.d0)
          !rho = 1.53d0
          q_loc(i,j,1) = rho
          q_loc(i,j,2) = vx*rho
          q_loc(i,j,3) = vy*rho
     enddo
   enddo

end subroutine
