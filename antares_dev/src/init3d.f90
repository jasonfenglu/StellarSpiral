subroutine init3d(q_loc)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR)
double precision::r_loc
integer::i,j,k
double precision::den,vx,vy,vz,ene,bxl,byl,bzl,bxr,byr,bzr

! implement your initial condition here
do k=1,ncell_loc(3)
  do j=1,ncell_loc(2)
    do i=1,ncell_loc(1)
      r_loc=dsqrt(x_loc(i)**2+y_loc(j)**2+z_loc(k)**2)
      den = 1.d0*dexp(-r_loc**2/2.d0)
      vx = 0.d0
      vy = 0.d0
      vz = 0.d0
      bxl = 0.75d0
      byl = 1.d0
      bzl = 0.d0
      bxr = 0.75d0
      byr = 1.d0
      bzr = 0.d0
     ene = 5.d0 
#ifdef MHD
     q_loc(i,j,k,1) = den   ! density
     q_loc(i,j,k,2) = den*vx  ! x momentum density
     q_loc(i,j,k,3) = den*vy  ! y momentum density
     q_loc(i,j,k,4) = den*vz  ! z momentum density
     q_loc(i,j,k,5) = bxl
     q_loc(i,j,k,6) = byl
     q_loc(i,j,k,7) = bzl
     q_loc(i,j,k,8) = ene
     q_loc(i,j,k,9) = bxr
     q_loc(i,j,k,10)= byr
     q_loc(i,j,k,11)= bzr
#else
     q_loc(i,j,1) = den   ! density
     q_loc(i,j,2) = den*vx  ! x momentum density
     q_loc(i,j,3) = den*vy  ! y momentum density
     q_loc(i,j,4) = 0.d0
#endif
    enddo
  enddo
enddo


end subroutine
