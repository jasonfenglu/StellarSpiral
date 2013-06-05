subroutine dt_all(dt_loc)
use common_params
implicit none
include 'mpif.h'
double precision::dt_loc,dtmin
integer::istatus(MPI_STATUS_SIZE), itag1,ierr,ireqs
double precision, dimension(:), allocatable:: r0dt,sdt,rdt,s0dt
integer::i

if(nprocs .gt. 1) then

allocate(r0dt(1:nprocs))
allocate(sdt(1))
allocate(rdt(1))
allocate(s0dt(1))

itag1=1
sdt(1)=dt_loc
!  send local dt to precessor 0 
if (myid.ne.0) then
       call MPI_Send(sdt,1,MPI_DOUBLE_PRECISION,0,itag1, &
             MPI_COMM_WORLD,ierr)
endif
!  receive dt from other processors 
if (myid.eq.0) then
   r0dt(1)=dt_loc
   do i=2,nprocs
       call MPI_Recv(r0dt(i),1,MPI_DOUBLE_PRECISION,i-1,itag1, &
             MPI_COMM_WORLD,istatus,ierr)
   enddo
!   determine global  dt 
   dtmin = r0dt(1)
!             write(6,*) 'dtmin=',dtmin, ' myid=',myid
   do i=2,nprocs
      dtmin= dmin1(dtmin,r0dt(i))
   enddo
   dt = dtmin
!   send dt to all processors 
!              write(6,*) 'Here dt=',dt, ' myid=',myid
   s0dt(1) = dt
   itag1   = 2
   do i=2,nprocs
     call MPI_Send(s0dt,1,MPI_DOUBLE_PRECISION,i-1,itag1, &
           MPI_COMM_WORLD,ierr)
   enddo
endif

if (myid.ne.0) then
  itag1 = 2
  call MPI_Recv(rdt,1,MPI_DOUBLE_PRECISION,0,itag1, &
           MPI_COMM_WORLD,istatus,ierr)
  dt = rdt(1)
endif

deallocate(r0dt)
deallocate(sdt)
deallocate(rdt)
deallocate(s0dt)

else
dt = dt_loc
endif
end subroutine
