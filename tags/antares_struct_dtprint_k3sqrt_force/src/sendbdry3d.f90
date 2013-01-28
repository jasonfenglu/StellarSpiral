subroutine sendbdry3d(q_loc)
use common_params
implicit none
include 'mpif.h'
double precision:: q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR)
integer:: i,j,k
double precision,dimension(:,:,:,:),allocatable::receive_buffer,send_buffer
integer::data_len
integer::istatus(MPI_STATUS_SIZE)
integer::ireqr,ireqs,ierr,itag1

allocate(send_buffer   (1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1:kbuf,NVAR))
allocate(receive_buffer(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1:kbuf,NVAR))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! send data to the right
itag1=3
data_len = (ncell_loc(1)+2*ibuf)*(ncell_loc(2)+2*jbuf)*kbuf*nvar
if (myid .ne. nprocs-1) then
  do k=1,kbuf 
    send_buffer(:,:,k,:)=q_loc(:,:,ncell_loc(3)-kbuf+k,:)
  enddo
  
  call MPI_ISend(send_buffer,data_len,MPI_DOUBLE_PRECISION,myid+1,itag1, &
               MPI_Comm_world,ireqs,ierr)
  call MPI_Wait(ireqs,istatus,ierr)
endif
! receive data from the left
if (myid .ne. 0) then
  call MPI_IRecv(receive_buffer,data_len,MPI_DOUBLE_PRECISION,myid-1,itag1, &
               MPI_COMM_WORLD,ireqr,ierr)
  call MPI_Wait(ireqr,istatus,ierr)
  
  do k=1,kbuf
    q_loc(:,:,k-kbuf,:)=receive_buffer(:,:,k,:)
  enddo
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! send data to the left 
itag1=4
if (myid .ne. 0) then
  do k=1,kbuf 
    send_buffer(:,:,k,:)=q_loc(:,:,k,:)
  enddo
  
  call MPI_ISend(send_buffer,data_len,MPI_DOUBLE_PRECISION,myid-1,itag1, &
               MPI_Comm_world,ireqs,ierr)
  call MPI_Wait(ireqs,istatus,ierr)
endif
! receive data from the right
if (myid .ne. nprocs-1) then
  call MPI_IRecv(receive_buffer,data_len,MPI_DOUBLE_PRECISION,myid+1,itag1, &
               MPI_COMM_WORLD,ireqr,ierr)
  call MPI_Wait(ireqr,istatus,ierr)

  do k=1,kbuf
    q_loc(:,:,ncell_loc(3)+k,:)=receive_buffer(:,:,k,:)
  enddo
endif

deallocate(send_buffer)
deallocate(receive_buffer)

end subroutine
