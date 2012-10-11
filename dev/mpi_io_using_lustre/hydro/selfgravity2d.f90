subroutine selfgravity2d(den)
#if NDIM==2
use common_params
implicit none
include "mpif.h"
include 'fftw_f77.i'

double precision::den(1:ncell_loc(1),1:ncell_loc(2))

double complex,dimension(:,:),allocatable::den_cmplx,den_work
double complex,dimension(:,:),allocatable::sgx_cmplx,sgy_cmplx
double precision,dimension(:,:),allocatable::send_buf1,send_buf2
double precision,dimension(:,:),allocatable::recv_buf1,recv_buf2
integer:: send_id1,send_id2,recv_id1,recv_id2,data_len
integer:: itag1,itag2,ireqs,ireqr,istatus,ierr
integer(8)::plan


allocate(den_cmplx(1:2*ncell_loc(1),1:2*ncell_loc(2)))
allocate( den_work(1:2*ncell_loc(1),1:2*ncell_loc(2)))
allocate(sgx_cmplx(1:2*ncell_loc(1),1:2*ncell_loc(2)))
allocate(sgy_cmplx(1:2*ncell_loc(1),1:2*ncell_loc(2)))
allocate(send_buf1(1:ncell_loc(1),1:ncell_loc(2)))
allocate(send_buf2(1:ncell_loc(1),1:ncell_loc(2)))
allocate(recv_buf1(1:ncell_loc(1),1:ncell_loc(2)))
allocate(recv_buf2(1:ncell_loc(1),1:ncell_loc(2)))

!==========   exchange data ================
if (nprocs .gt. 1) then
itag1=1
send_id1 = int(floor(float(myid)/2.))
recv_id1 = 2*myid
recv_id2 = 2*myid+1
data_len = ncell_loc(1)*ncell_loc(2)
den_cmplx = dcmplx(0.d0,0.d0)

send_buf1 = den(1:ncell_loc(1),1:ncell_loc(2))


call MPI_ISend(send_buf1,data_len,MPI_DOUBLE_PRECISION,send_id1,itag1, &
               MPI_COMM_WORLD,ireqs,ierr)

if(myid .lt. nprocs/2) then
  itag1=1
  call MPI_IRecv(recv_buf1,data_len,MPI_DOUBLE_PRECISION,recv_id1,itag1, &
               MPI_COMM_WORLD,ireqr,ierr)
  call MPI_Wait(ireqr,istatus,ierr)
  itag1=1
  call MPI_IRecv(recv_buf2,data_len,MPI_DOUBLE_PRECISION,recv_id2,itag1, &
               MPI_COMM_WORLD,ireqr,ierr)
  call MPI_Wait(ireqr,istatus,ierr)
  den_cmplx(1:ncell_loc(1),1:ncell_loc(2))=dcmplx(recv_buf1,0.d0)
  den_cmplx(1:ncell_loc(1),ncell_loc(2)+1:2*ncell_loc(2))=dcmplx(recv_buf2,0.d0)
endif
call MPI_Wait(ireqs,istatus,ierr)
else
den_cmplx(1:ncell_loc(1),1:ncell_loc(2)) = dcmplx(den(1:ncell_loc(1),1:ncell_loc(2)),0.d0)
endif  ! end if nprocs

if(myid .eq. 0) then
 write(*,*) "calculating selfgravity"
endif
!========== do FFT to density =======================================
call fftw2d_f77_mpi_create_plan(plan,MPI_COMM_WORLD,2*ncell(1),2*ncell(2),FFTW_FORWARD,FFTW_ESTIMATE)
call fftwnd_f77_mpi(plan,1,den_cmplx,den_work,1,FFTW_NORMAL_ORDER)
call fftwnd_f77_mpi_destroy_plan(plan)
!=================================================
sgx_cmplx = sgxker*den_cmplx
sgy_cmplx = sgyker*den_cmplx

call fftw2d_f77_mpi_create_plan(plan,MPI_COMM_WORLD,2*ncell(1),2*ncell(2),FFTW_BACKWARD,FFTW_ESTIMATE)
call fftwnd_f77_mpi(plan,1,sgx_cmplx,den_work,1,FFTW_NORMAL_ORDER)
call fftwnd_f77_mpi_destroy_plan(plan)

call fftw2d_f77_mpi_create_plan(plan,MPI_COMM_WORLD,2*ncell(1),2*ncell(2),FFTW_BACKWARD,FFTW_ESTIMATE)
call fftwnd_f77_mpi(plan,1,sgy_cmplx,den_work,1,FFTW_NORMAL_ORDER)
call fftwnd_f77_mpi_destroy_plan(plan)

!========== send x-force back to correct processors =====
if(nprocs .gt. 1) then
recv_id1 = int(floor(float(myid+nprocs)/2.))
send_id1 = 2*myid-nprocs
send_id2 = 2*myid-nprocs+1

if(myid .ge. nprocs/2) then
  send_buf1 = dreal(sgx_cmplx(ncell_loc(1)+1:2*ncell_loc(1),1:ncell_loc(2)))/dble(4*ncell(1)*ncell(2))*GravConst
  send_buf2 = dreal(sgx_cmplx(ncell_loc(1)+1:2*ncell_loc(1),ncell_loc(2)+1:2*ncell_loc(2)))/dble(4*ncell(1)*ncell(2))*GravConst
  itag1=1
  call MPI_ISend(send_buf1,data_len,MPI_DOUBLE_PRECISION,send_id1,itag1, &
                MPI_Comm_world,ireqs,ierr)
  call MPI_Wait(ireqs,istatus,ierr)
  itag2=2
  call MPI_ISend(send_buf2,data_len,MPI_DOUBLE_PRECISION,send_id2,itag2, &
               MPI_Comm_world,ireqs,ierr)
  call MPI_Wait(ireqs,istatus,ierr)
endif

if(mod(myid,2) .eq. 0) then
  itag1=1
  call MPI_IRecv(sgx,data_len,MPI_DOUBLE_PRECISION,recv_id1,itag1, &
               MPI_COMM_WORLD,ireqr,ierr)
  call MPI_Wait(ireqr,istatus,ierr)
endif

if(mod(myid,2) .eq. 1) then
  itag2=2
  call MPI_IRecv(sgx,data_len,MPI_DOUBLE_PRECISION,recv_id1,itag2, &
               MPI_COMM_WORLD,ireqr,ierr)
  call MPI_Wait(ireqr,istatus,ierr)
endif
call MPI_barrier(MPI_COMM_WORLD,ierr)

else
sgx=dreal(sgx_cmplx(ncell_loc(1)+1:2*ncell_loc(1),ncell_loc(2)+1:2*ncell_loc(2)))/dble(4*ncell(1)*ncell(2))*GravConst
endif  ! end if nprocs

!========== send y-force back to correct processors =====
if(nprocs .gt. 1) then
recv_id1 = int(floor(float(myid+nprocs)/2.))
send_id1 = 2*myid-nprocs
send_id2 = 2*myid-nprocs+1

if(myid .ge. nprocs/2) then
  send_buf1 =dreal(sgy_cmplx(ncell_loc(1)+1:2*ncell_loc(1),1:ncell_loc(2)))/dble(4*ncell(1)*ncell(2))*GravConst
  send_buf2 =dreal(sgy_cmplx(ncell_loc(1)+1:2*ncell_loc(1),ncell_loc(2)+1:2*ncell_loc(2)))/dble(4*ncell(1)*ncell(2))*GravConst
  itag1=1
  call MPI_ISend(send_buf1,data_len,MPI_DOUBLE_PRECISION,send_id1,itag1, &
                MPI_Comm_world,ireqs,ierr)
  call MPI_Wait(ireqs,istatus,ierr)
  itag2=2
  call MPI_ISend(send_buf2,data_len,MPI_DOUBLE_PRECISION,send_id2,itag2, &
               MPI_Comm_world,ireqs,ierr)
  call MPI_Wait(ireqs,istatus,ierr)
endif

if(mod(myid,2) .eq. 0) then
  itag1=1
  call MPI_IRecv(sgy,data_len,MPI_DOUBLE_PRECISION,recv_id1,itag1, &
               MPI_COMM_WORLD,ireqr,ierr)
  call MPI_Wait(ireqr,istatus,ierr)
endif

if(mod(myid,2) .eq. 1) then
  itag2=2
  call MPI_IRecv(sgy,data_len,MPI_DOUBLE_PRECISION,recv_id1,itag2, &
               MPI_COMM_WORLD,ireqr,ierr)
  call MPI_Wait(ireqr,istatus,ierr)
endif
call MPI_barrier(MPI_COMM_WORLD,ierr)

else
sgy=dreal(sgy_cmplx(ncell_loc(1)+1:2*ncell_loc(1),ncell_loc(2)+1:2*ncell_loc(2)))/dble(4*ncell(1)*ncell(2))*GravConst
endif  ! end if nprocs

deallocate(den_cmplx)
deallocate( den_work)
deallocate(sgx_cmplx)
deallocate(sgy_cmplx)
deallocate(send_buf1)
deallocate(send_buf2)
deallocate(recv_buf1)
deallocate(recv_buf2)
#endif
end subroutine
