subroutine sgkernel3d()
#if NDIM==3
use common_params
implicit none
include 'mpif.h'
include 'fftw_f77.i'
double complex,dimension(:,:,:),allocatable::sgker_work
double precision,dimension(:),allocatable::x_ker_loc,y_ker_loc,z_ker_loc
double precision::r_ker_loc,hdx,hdy,hdz,xl,yl,zl,dvolume
integer:: i,j,k
integer(8)::plan

allocate(x_ker_loc(1:2*ncell_loc(1)))
allocate(y_ker_loc(1:2*ncell_loc(2)))
allocate(z_ker_loc(1:2*ncell_loc(3)))
allocate(sgker_work(1:2*ncell_loc(1),1:2*ncell_loc(2),1:2*ncell_loc(3)))

do i=1,2*ncell_loc(1)
  x_ker_loc(i)=2*xrange(1)+(dble(i)-0.5d0)*dx 
enddo

do i=1,2*ncell_loc(2)
  y_ker_loc(i)=2*yrange(1)+(dble(i)-0.5d0)*dy
enddo

do i=1,2*ncell_loc(3)
  z_ker_loc(i)=2*zrange(1)+(dble(i)-0.5d0+dble(myid*2*ncell_loc(3)))*dz
enddo

hdx=0.5d0*(x_loc(2)-x_loc(1))
hdy=0.5d0*(y_loc(2)-y_loc(1))
hdz=0.5d0*(z_loc(2)-z_loc(1))
dvolume = dx*dy*dz

sgxker=dcmplx(0.d0,0.d0)
sgyker=dcmplx(0.d0,0.d0)
sgzker=dcmplx(0.d0,0.d0)

if(myid .eq. 0) then
   write(*,*) "Filling 3D kernal matrix for self-gravity......"
endif

do k=1,2*ncell_loc(3)
 do j=1,2*ncell_loc(2)
  do i=1,2*ncell_loc(1)
    zl = z_ker_loc(k)-hdz
    yl = y_ker_loc(j)-hdy
    xl = x_ker_loc(i)-hdx
    
    r_ker_loc = dsqrt(yl**2+xl**2+zl**2)
    sgxker(i,j,k)=dcmplx(-xl/r_ker_loc**3*dvolume,0.d0)
    sgyker(i,j,k)=dcmplx(-yl/r_ker_loc**3*dvolume,0.d0)
    sgzker(i,j,k)=dcmplx(-zl/r_ker_loc**3*dvolume,0.d0)

      if (r_ker_loc .lt. 0.5d0*dx) then
        sgxker(i,j,k) = dcmplx(0.d0,0.d0)
        sgyker(i,j,k) = dcmplx(0.d0,0.d0)
        sgzker(i,j,k) = dcmplx(0.d0,0.d0)
      endif
!!! The following is obtained from integrating over a single cell,
!!! in principle, this should be more accurate
!!!!!!! preparing Kx
!       xtemp1 = (yr-hdx)/dabs(xr-hdx)
!       xtemp2 = (yr+hdx)/dabs(xr-hdx)
!       xtemp3 = (yr+hdx)/dabs(xr+hdx)
!       xtemp4 = (yr-hdx)/dabs(xr+hdx)
!      sgxker(i,j) = dcmplx(asinh(-xtemp1)-asinh(-xtemp2)-asinh(-xtemp4)+asinh(-xtemp3),0.d0)
!!!!!! preparing Ky
!       ytemp1 = (xr-hdy)/dabs(yr-hdy)
!       ytemp2 = (xr+hdy)/dabs(yr-hdy)
!       ytemp3 = (xr+hdy)/dabs(yr+hdy)
!       ytemp4 = (xr-hdy)/dabs(yr+hdy)
!      sgyker(i,j) =dcmplx(asinh(-ytemp1)-asinh(-ytemp2)-asinh(-ytemp4)+asinh(-ytemp3),0.d0)
    enddo
  enddo
enddo
!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!
  if (myid .eq. 0) then
     write(*,*) "Complete!"
     write(*,*) "sgkernel2d.f90:doing 3D-FFT to kernal"
  endif

call fftw3d_f77_mpi_create_plan(plan,MPI_COMM_WORLD,2*ncell(1),2*ncell(2),2*ncell(3), FFTW_FORWARD,FFTW_ESTIMATE)
call fftwnd_f77_mpi(plan,1,sgxker,sgker_work,1,FFTW_NORMAL_ORDER)
call fftwnd_f77_mpi_destroy_plan(plan)

call fftw3d_f77_mpi_create_plan(plan,MPI_COMM_WORLD,2*ncell(1),2*ncell(2),2*ncell(3), FFTW_FORWARD,FFTW_ESTIMATE)
call fftwnd_f77_mpi(plan,1,sgyker,sgker_work,1,FFTW_NORMAL_ORDER)
call fftwnd_f77_mpi_destroy_plan(plan)

call fftw3d_f77_mpi_create_plan(plan,MPI_COMM_WORLD,2*ncell(1),2*ncell(2),2*ncell(3), FFTW_FORWARD,FFTW_ESTIMATE)
call fftwnd_f77_mpi(plan,1,sgzker,sgker_work,1,FFTW_NORMAL_ORDER)
call fftwnd_f77_mpi_destroy_plan(plan)

  if (myid .eq. 0) then
     write(*,*) "sgkernel3d.f90:3D kernel FFT Done!!"
  endif

deallocate(x_ker_loc)
deallocate(y_ker_loc)
deallocate(z_ker_loc)
deallocate(sgker_work)
#endif
end subroutine
