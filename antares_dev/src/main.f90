program antares
use common_params
use simcontroll
use STELLARDISK,only:INIT_STELLARDISK,k3sqrtlog
implicit none
include 'mpif.h'
include 'fftw_f77.i'

character(len=20)::infile='patch/para.nml' ! filename of namelist
character(len=8)::flnm
character(len=20)::dsetname
integer::ierr ! error flag for MPI
double precision::dt_loc,toutput ! local time step, next point for data output
integer::fnum=0,fstart! number of file number
integer::i,j,k,step

#if NDIM==2
double precision,dimension(:,:,:),allocatable::q_loc,temp1_loc,temp2_loc
double precision,dimension(:,:),allocatable::den_temp
#ifdef MHD
double precision,dimension(:,:,:),allocatable::temp3_loc
#endif

#ifdef FORCE
double precision,dimension(:,:),allocatable::fx,fy  ! external force
#endif

#ifdef GRAVITY
double precision,dimension(:,:),allocatable::den 
#endif
#endif

#if NDIM==3
double precision,dimension(:,:,:,:),allocatable::q_loc,temp1_loc,temp2_loc

#ifdef FORCE
double precision,dimension(:,:,:),allocatable::fx,fy,fz
#endif

#ifdef GRAVITY
double precision,dimension(:,:,:),allocatable::den 
#endif
#endif




namelist/run_params/ncell,xrange,yrange,zrange, &
                    tend,CFL,&
                    dtout,snd,gam,GravConst,&
                    fstart

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

#ifdef VERBOSE
write(*,'("node",I2,"/",I2,"   ready!!")') myid+1,nprocs
#endif

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
open(1,file=infile)
read(1,NML=run_params)
close(1)

!!init by ccfeng
call init_simcon(dtout,tend)
call INIT_STELLARDISK(ncell(2),xrange(2))
if(myid .eq. 0)CALL k3sqrtlog


#ifdef VERBOSE
if(myid .eq. 0) then
  write(*,*) "NDIM: ",NDIM
  write(*,*) "ncell: ",ncell
  write(*,*) "xrange: ",xrange
  write(*,*) "yrange: ",yrange
  write(*,*) "zrange: ",zrange
  write(*,*) "tend: ",tend
  write(*,*) "CFL: ",CFL
  write(*,*) "dtout: ",dtout
  write(*,*) "snd: ",snd
  write(*,*) "gam: ",gam
  write(*,*) "fstart: ",fstart
#ifdef GRAVITY
  write(*,*) "GravConst: ",GravConst
#endif
endif
#endif
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#if NDIM==2

ncell_loc(1) = ncell(1)
ncell_loc(2) = ncell(2)/nprocs
ncell_loc(3) = 0
allocate(x_loc(1-ibuf:ncell_loc(1)+ibuf))
allocate(y_loc(1-jbuf:ncell_loc(2)+jbuf))
allocate(x(1-ibuf:ncell(1)+ibuf))
allocate(y(1-jbuf:ncell(2)+jbuf))
allocate(q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR))
allocate(temp1_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR))
allocate(temp2_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR))
allocate(den_temp(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))

#ifdef FORCE
allocate(fx(1:ncell_loc(1),1:ncell_loc(2)))
allocate(fy(1:ncell_loc(1),1:ncell_loc(2)))
#endif

#ifdef GRAVITY
allocate(den(1:ncell_loc(1),1:ncell_loc(2)))
allocate(sgx(1:ncell_loc(1),1:ncell_loc(2)))
allocate(sgy(1:ncell_loc(1),1:ncell_loc(2)))
allocate(sgxker(1:2*ncell_loc(1),1:2*ncell_loc(2)))
allocate(sgyker(1:2*ncell_loc(1),1:2*ncell_loc(2)))
#endif

dx = (xrange(2)-xrange(1))/dble(ncell(1))
dy = (yrange(2)-yrange(1))/dble(ncell(2))

do i = 1-ibuf, ncell_loc(1)+ibuf
    x_loc(i) = xrange(1)+(dble(i)-0.5d0)*dx
enddo

do i = 1-jbuf, ncell_loc(2)+jbuf
    y_loc(i) = yrange(1)+(dble(i)-0.5d0+dble(myid*ncell_loc(2)))*dy
enddo

do i=1-ibuf, ncell(1)+ibuf
   x(i) = xrange(1)+(dble(i)-0.5d0)*dx
enddo

do i=1-jbuf, ncell(2)+jbuf
   y(i) = yrange(1)+(dble(i)-0.5d0)*dy
enddo

if(fstart .eq. 0) then
  t = 0.d0
  step = 0
  call init2d(q_loc)  !! fresh start
  call bdry2d(q_loc)
  if(nprocs .gt. 1) then
    call sendbdry2d(q_loc)
  endif
  call write_collective2d(q_loc,fnum)
#ifdef GRAVITY
  call sgkernel2d() !!prepare the complex kernel FFT
  if(myid .eq. 0) then
  write(*,*) "calculating self-gravity for the initial condition......"
  endif
  den=q_loc(1:ncell_loc(1),1:ncell_loc(2),1)
  call selfgravity2d(den)  !!real force of selfgravity
  if(myid .eq. 0) then
    write(*,*) "done"
    write(*,*) "writing initial selfgravity into files......"
  endif
  call write_gravity2d(fnum)
  if(myid .eq. 0) then
    write(*,*) "done."
    write(*,*) "start main loop......"
  endif
#endif
else   ! not fresh start
  write(*,*) "loading ...."
  call read_collective2d(q_loc,fstart) !!read data from fstart
  if(nprocs .gt. 1) then
     call sendbdry2d(q_loc)
  endif
#ifdef GRAVITY
  call sgkernel2d() !!prepare the complex kernel FFT
#endif
  fnum = fstart
endif

call dt_loc2d(q_loc,dt_loc)
call dt_all(dt_loc)
toutput=t+dtout

!!!! main loop
do while(t .lt. tend)
t=t+dt
step = step + 1
if(myid .eq. 0) then
   write(*,*) "(t,dt)=",step,t,dt
endif

#ifdef RK2
! RK1
call riemann2d(q_loc,q_loc,temp1_loc,1)  ! 1: x
call riemann2d(q_loc,temp1_loc,temp1_loc,2) ! 2: y


#ifdef FORCE
  call force2d(q_loc,fx,fy)  ! external force
  do j=1,ncell_loc(2)
    do i=1,ncell_loc(1)
      temp1_loc(i,j,2) = temp1_loc(i,j,2)+fx(i,j)*dt*q_loc(i,j,1)
      temp1_loc(i,j,3) = temp1_loc(i,j,3)+fy(i,j)*dt*q_loc(i,j,1)
!!!! For non-isothermal simulation, one needs to take into account the energy equation 
    enddo
  enddo
#endif


call bdry2d(temp1_loc)
if(nprocs .gt. 1) then
   call sendbdry2d(temp1_loc)
endif


! RK2
call riemann2d(temp1_loc,temp1_loc,temp2_loc,1)
call riemann2d(temp1_loc,temp2_loc,temp2_loc,2)

#ifdef FORCE
  call force2d(temp1_loc,fx,fy)
  do j=1,ncell_loc(2)
    do i=1,ncell_loc(1)
      temp2_loc(i,j,2)=temp2_loc(i,j,2)+fx(i,j)*dt*temp1_loc(i,j,1)
      temp2_loc(i,j,3)=temp2_loc(i,j,3)+fy(i,j)*dt*temp1_loc(i,j,1)
!!!! For non-isothermal simulation, one needs to take into account the energy equation 
    enddo
  enddo 
#endif

call bdry2d(temp2_loc)
if(nprocs .gt. 1) then
   call sendbdry2d(temp2_loc)
endif

   q_loc(1:ncell_loc(1),1:ncell_loc(2),:)=0.5d0*(q_loc(1:ncell_loc(1),1:ncell_loc(2),:)+temp2_loc(1:ncell_loc(1),1:ncell_loc(2),:))

#endif
!!!!! end of RK2

#ifdef CTU
  den_temp = q_loc(:,:,1)
  call force2d(q_loc,fx,fy)
  call riemann2dCTU(q_loc)
  !call write_collective2d(q_loc,step)
#ifdef FORCE
  do j=1,ncell_loc(2)
    do i=1,ncell_loc(1)
      q_loc(i,j,2)=q_loc(i,j,2)+fx(i,j)*dt*den_temp(i,j)
      q_loc(i,j,3)=q_loc(i,j,3)+fy(i,j)*dt*den_temp(i,j)
    enddo
  enddo
#endif
#endif

call bdry2d(q_loc)
if(nprocs .gt. 1) then
   call sendbdry2d(q_loc)
endif


if (t .ge. toutput) then
  fnum = fnum+1
  call write_collective2d(q_loc,fnum)
#ifdef GRAVITY
  call write_gravity2d(fnum)
#endif
  if(myid .eq. 0) then
     write(*,*) myid,"output data:",fnum
  endif
  if(toutput+dtout .gt. tend) then
    toutput = tend
  else
    toutput = toutput + dtout
  endif
endif


call dt_loc2d(q_loc,dt_loc)
call dt_all(dt_loc)
!!add by ccfeng
if(dt.lt.1d-10)then
        if(myid.eq.0)print *,'dt too small, quit'
        call write_collective2d(q_loc,fnum+1)
        call MPI_BARRIER(MPI_COMM_WORLD)
        exit
endif

if(t+dt .ge. toutput) then
   dt = toutput-t
endif

if(t+dt .ge. tend) then
   dt = tend-t
endif

enddo  ! end of main loop

deallocate(x_loc)
deallocate(y_loc)
deallocate(x)
deallocate(y)
deallocate(q_loc)
deallocate(temp1_loc)
deallocate(temp2_loc)

#ifdef FORCE
deallocate(fx)
deallocate(fy)
#endif

#ifdef GRAVITY
deallocate(den)
deallocate(sgx)
deallocate(sgy)
deallocate(sgxker)
deallocate(sgyker)
#endif
#endif

!!!!!!!!!!!!!!!!!! 3D solver
#if NDIM==3
ncell_loc(1)=ncell(1)
ncell_loc(2)=ncell(2)
ncell_loc(3)=ncell(3)/nprocs

allocate(x_loc(1-ibuf:ncell_loc(1)+ibuf))
allocate(y_loc(1-jbuf:ncell_loc(2)+jbuf))
allocate(z_loc(1-kbuf:ncell_loc(3)+kbuf))
allocate(x(1-ibuf:ncell(1)+ibuf))
allocate(y(1-jbuf:ncell(2)+jbuf))
allocate(z(1-kbuf:ncell(3)+kbuf))
allocate(q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR))
allocate(temp1_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR))
allocate(temp2_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR))

#ifdef FORCE
allocate(fx(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3)))
allocate(fy(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3)))
allocate(fz(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3)))
#endif

#ifdef GRAVITY
allocate(den(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3)))
allocate(sgx(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3)))
allocate(sgy(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3)))
allocate(sgz(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3)))
allocate(sgxker(1:2*ncell_loc(1),1:2*ncell_loc(2),1:2*ncell_loc(3)))
allocate(sgyker(1:2*ncell_loc(1),1:2*ncell_loc(2),1:2*ncell_loc(3)))
allocate(sgzker(1:2*ncell_loc(1),1:2*ncell_loc(2),1:2*ncell_loc(3)))
#endif

dx=(xrange(2)-xrange(1))/dble(ncell(1))
dy=(yrange(2)-yrange(1))/dble(ncell(2))
dz=(zrange(2)-zrange(1))/dble(ncell(3))

do i=1-ibuf,ncell_loc(1)+ibuf
  x_loc(i)=xrange(1)+(dble(i)-0.5d0)*dx
enddo

do j=1-jbuf,ncell_loc(2)+jbuf
  y_loc(j)=yrange(1)+(dble(j)-0.5d0)*dy
enddo

do k=1-kbuf,ncell_loc(3)+kbuf
  z_loc(k)=zrange(1)+(dble(k)-0.5d0+dble(myid*ncell_loc(3)))*dz
enddo

do i=1-ibuf,ncell_loc(1)+ibuf
  x(i)=xrange(1)+(dble(i)-0.5d0)*dx
enddo

do j=1-jbuf,ncell_loc(2)+jbuf
  y(j)=yrange(1)+(dble(j)-0.5d0)*dy
enddo

do k=1-kbuf,ncell_loc(3)+kbuf
  z(k)=zrange(1)+(dble(k)-0.5d0)*dz
enddo

if(fstart .eq. 0)then
  t=0.d0
  call init3d(q_loc)
  call bdry3d(q_loc)
  if(nprocs .gt. 1) then
    call sendbdry3d(q_loc)
  endif
  call write_collective3d(q_loc,fnum)
#ifdef GRAVITY
  call sgkernel3d() !! prepare the complex kernel FFT 
  if(myid .eq. 0) then
    write(*,*) "calculating self-gravity for the initial condition...."
  endif
  den=q_loc(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3),1)
  call selfgravity3d(den)
  if(myid .eq. 0) then
    write(*,*) "done."
    write(*,*) "writing initial selfgravity into files......"
  endif
  call write_gravity3d(fnum)
  if(myid .eq. 0) then
    write(*,*) "done."
    write(*,*) "start main loop......"
  endif
#endif
else
endif

call dt_loc3d(q_loc,dt_loc)
call dt_all(dt_loc)
toutput=t+dtout


!!! main loop
do while(t .lt. tend)
t=t+dt
if(myid .eq. 0) then
  write(*,*) "(t,dt)",t,dt
endif
! RK1
call riemann3d(q_loc,q_loc    ,temp1_loc,1)
call riemann3d(q_loc,temp1_loc,temp1_loc,2)
call riemann3d(q_loc,temp1_loc,temp1_loc,3)

#ifdef FORCE
   call force3d(q_loc,fx,fy,fz)
   do k=1,ncell_loc(3)
     do j=1,ncell_loc(2)
       do i=1,ncell_loc(1)
          temp1_loc(i,j,k,2)=temp1_loc(i,j,k,2)+fx(i,j,k)*dt*q_loc(i,j,k,1)
          temp1_loc(i,j,k,3)=temp1_loc(i,j,k,3)+fy(i,j,k)*dt*q_loc(i,j,k,1)
          temp1_loc(i,j,k,4)=temp1_loc(i,j,k,4)+fz(i,j,k)*dt*q_loc(i,j,k,1)
       enddo
     enddo
  enddo
#endif

call bdry3d(temp1_loc)
if(nprocs .gt. 1) then
   call sendbdry3d(temp1_loc)
endif

!! RK2
call riemann3d(temp1_loc,temp1_loc,temp2_loc,1)
call riemann3d(temp1_loc,temp2_loc,temp2_loc,2)
call riemann3d(temp1_loc,temp2_loc,temp2_loc,3)

#ifdef FORCE
   call force3d(temp1_loc,fx,fy,fz)
   do k=1,ncell_loc(3)
     do j=1,ncell_loc(2)
       do i=1,ncell_loc(1)
          temp2_loc(i,j,k,2)=temp2_loc(i,j,k,2)+fx(i,j,k)*dt*temp1_loc(i,j,k,1)
          temp2_loc(i,j,k,3)=temp2_loc(i,j,k,3)+fy(i,j,k)*dt*temp1_loc(i,j,k,1)
          temp2_loc(i,j,k,4)=temp2_loc(i,j,k,4)+fz(i,j,k)*dt*temp1_loc(i,j,k,1)
       enddo
     enddo
   enddo
#endif


call bdry3d(temp2_loc)
if(nprocs .gt. 1) then
   call sendbdry3d(temp2_loc)
endif

q_loc(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3),:)=0.5d0*(q_loc(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3),:)+temp2_loc(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3),:))

call bdry3d(q_loc)
if(nprocs .gt. 1) then
   call sendbdry3d(q_loc)
endif

if (t .ge. toutput) then
  fnum = fnum+1
  call write_collective3d(q_loc,fnum)
#ifdef GRAVITY
  call write_gravity3d(fnum)
#endif
  if(myid .eq. 0) then
     write(*,*) myid,"output data:",fnum
  endif
  if(toutput+dtout .gt. tend) then
    toutput = tend
  else
    toutput = toutput + dtout
 endif
endif


call dt_loc3d(q_loc,dt_loc)
call dt_all(dt_loc)

if(t+dt .ge. toutput) then
   dt = toutput-t
endif

if(t+dt .ge. tend) then
   dt = tend-t
endif

enddo

deallocate(x_loc)
deallocate(y_loc)
deallocate(z_loc)
deallocate(x)
deallocate(y)
deallocate(z)
deallocate(q_loc)
deallocate(temp1_loc)
deallocate(temp2_loc)

#ifdef FORCE
deallocate(fx)
deallocate(fy)
deallocate(fz)
#endif

#ifdef GRAVITY
deallocate(den)
deallocate(sgx)
deallocate(sgy)
deallocate(sgz)
deallocate(sgxker)
deallocate(sgyker)
deallocate(sgzker)
#endif
#endif

call MPI_FINALIZE(ierr)

end program
