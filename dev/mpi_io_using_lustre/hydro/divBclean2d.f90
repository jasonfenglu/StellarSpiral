subroutine divBclean2d(q_loc1,q_loc2)
use common_params
implicit none
double precision::q_loc1(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::q_loc2(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision, dimension(:,:),allocatable:: vB1,vB2,Ez
integer:: i,j

allocate(vB1(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))
allocate(vB2(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))
allocate( Ez(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))

vB1 = q_loc1(:,:,6)*q_loc1(:,:,2)/q_loc1(:,:,1)-q_loc1(:,:,5)*q_loc1(:,:,3)/q_loc1(:,:,1)
vB2 = q_loc2(:,:,6)*q_loc2(:,:,2)/q_loc2(:,:,1)-q_loc2(:,:,5)*q_loc2(:,:,3)/q_loc2(:,:,1)
Ez  = -0.5d0*(vB1+vB2)
!Ez = -vB2

do j=1,ncell_loc(2)
  do i=1,ncell_loc(1)
   q_loc2(i,j,5) = q_loc1(i,j,5)-(Ez(i,j+1)-Ez(i,j-1))*dt/(2.d0*dx)
   q_loc2(i,j,6) = q_loc1(i,j,6)+(Ez(i+1,j)-Ez(i-1,j))*dt/(2.d0*dx)
  enddo
enddo

deallocate(vB1)
deallocate(vB2)
deallocate(Ez)

end subroutine
