subroutine bdry2d(q_loc)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
integer::i,j,k


!!!! implement your boundary condition here

!!! default boundary condition: zero gradient extrapolation
!do j=1,ncell_loc(2)
  do i=1,ibuf
    !do k=1,NVAR
    !   q_loc(1-i,j,k) = q_loc(1,j,k)
    !   q_loc(ncell_loc(1)+i,j,k) = q_loc(ncell_loc(1),j,k)
       q_loc(1-i,:,:) = q_loc(1,:,:)
       q_loc(ncell_loc(1)+i,:,:) = q_loc(ncell_loc(1),:,:)
     
    !enddo
  enddo
!enddo


if(myid .eq. 0) then
!   do i=1,ncell_loc(1)
     do j=1,jbuf
 !      do k=1,NVAR
          q_loc(:,1-j,:) = q_loc(:,1,:)
 !      enddo
     enddo
!   enddo
endif

if(myid .eq. nprocs-1) then
!   do i=1,ncell_loc(1)
     do j=1,jbuf
 !      do k=1,NVAR
          q_loc(:,ncell_loc(2)+j,:)=q_loc(:,ncell_loc(2),:)
 !      enddo
     enddo
!   enddo
endif

end subroutine
