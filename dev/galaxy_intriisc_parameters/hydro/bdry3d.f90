subroutine bdry3d(q_loc)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR)
integer::i,j,k,l


!!!! implement your boundary condition here

!!! default boundary condition: zero gradient extrapolation
do i=1,ibuf
  !do k=1,ncell_loc(3)
  !  do j=1,ncell_loc(2)
  !    do l=1,NVAR
        q_loc(1-i,:,:,:) = q_loc(1,:,:,:)
        q_loc(ncell_loc(1)+i,:,:,:) = q_loc(ncell_loc(1),:,:,:)
  !    enddo
  !  enddo
  !enddo
enddo

do j=1,jbuf
  !do k=1,ncell_loc(3)
  !  do i=1,ncell_loc(1)
  !    do l=1,NVAR
        q_loc(:,1-j,:,:) = q_loc(:,1,:,:)
        q_loc(:,ncell_loc(2)+j,:,:) = q_loc(:,ncell_loc(2),:,:)
  !    enddo 
  !  enddo
  !enddo
enddo


if(myid .eq. 0) then
do k=1,kbuf
   !do j=1,ncell_loc(2)
   !  do i=1,ncell_loc(1)
   !    do l=1,NVAR
          q_loc(:,:,1-k,:) = q_loc(:,:,1,:)
   !    enddo
   !  enddo
   !enddo
enddo
endif

if(myid .eq. nprocs-1) then
do k=1,kbuf
   !do j=1,ncell_loc(2)
   !  do i=1,ncell_loc(1)
   !    do l=1,NVAR
          q_loc(:,:,ncell_loc(3)+k,:)=q_loc(:,:,ncell_loc(3),:)
   !    enddo
   !  enddo
   !enddo
enddo
endif

end subroutine
