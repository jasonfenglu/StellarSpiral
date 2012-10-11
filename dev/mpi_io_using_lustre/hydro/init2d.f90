subroutine init2d(q_loc)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::r_loc
integer::i,j
double precision::den,vx,vy,vz,ene,bx,by,bz

! implement your initial condition here

do j=1,ncell_loc(2)
  do i=1,ncell_loc(1)
   r_loc=dsqrt(x_loc(i)**2+y_loc(j)**2)
   den = 1.d0*dexp(-r_loc**2/2.d0)
    vx = 0.d0
    vy = 0.d0
    vz = 0.d0
    bx = 0.75d0
    by = 1.d0
    bz = 0.d0
   ene = 5.d0 
#ifdef MHD
  q_loc(i,j,1) = den   ! density
  q_loc(i,j,2) = den*vx  ! x momentum density
  q_loc(i,j,3) = den*vy  ! y momentum density
  q_loc(i,j,4) = den*vz  ! z momentum density
  q_loc(i,j,5) = bx
  q_loc(i,j,6) = by
  q_loc(i,j,7) = bz
  q_loc(i,j,8) = ene
#else
  q_loc(i,j,1) = den   ! density
  q_loc(i,j,2) = den*vx  ! x momentum density
  q_loc(i,j,3) = den*vy  ! y momentum density
  q_loc(i,j,4) = 0.d0
#endif
  enddo
enddo


end subroutine
