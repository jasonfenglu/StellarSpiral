subroutine output3d(q_loc,nx_global,ny_global,nz_global,nx_loc,ny_loc,nz_loc,inputbufx,inputbufy,inputbufz,bufx,bufy,bufz,flnm,dsetname,flag)
use hdf5
use common_params
implicit none
include 'mpif.h'

integer::nx_global,ny_global,nz_global,nx_loc,ny_loc,nz_loc,flag,bufx,inputbufx,bufy,inputbufy,bufz,inputbufz
double precision::q_loc(1-inputbufx:nx_loc+inputbufx,1-inputbufy:ny_loc+inputbufy,1-inputbufz:nz_loc+inputbufz)
double precision, dimension(:,:,:),allocatable::q_temp
character(len=20):: dsetname
character(len=8 ):: flnm  ! output filename

integer(HID_T)    ::file_id
integer(HID_T)    ::dset_id
integer(HID_T)    ::filespace
integer(HID_T)    ::mspace_id
integer(HID_T)    ::plist_id
integer(HSIZE_T),dimension(3) :: dimsf
integer(HSIZE_T),dimension(3) :: dimsf_loc
integer(HSIZE_T),dimension(3) :: istart, istride, icount, iblock
integer           ::error,ierr


!! create file
call h5open_f(error)
call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,error)
  if (flag .eq. 1) then
#ifdef VERBOSE
    if(myid .eq. 0) then
      write(*,*) "creating new file:",flnm
      write(*,*) "writing dataset ",dsetname, " into the file ", flnm
    endif
#endif
    call h5fcreate_f(flnm, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
  else
#ifdef VERBOSE
    if(myid .eq. 0) then
      write(*,*) "writing dataset ",dsetname, " into the file ", flnm
    endif
#endif
    call h5fopen_f(flnm, H5F_ACC_RDWR_F,file_id, error, access_prp=plist_id)
  endif
call h5pclose_f(plist_id, error)

!! prepare date to write 

dimsf(1) = nx_global+2*bufx
dimsf(2) = ny_global+2*bufy
dimsf(3) = nz_global+2*bufz

if(nz_global .ne. nz_loc) then
   dimsf_loc(1) = nx_loc+2*bufx
   dimsf_loc(2) = ny_loc+2*bufy
   dimsf_loc(3) = nz_loc

     istart(1) = 0 
     istart(2) = 0
     istart(3) = myid*nz_loc+bufz

   istride(1) = 1
   istride(2) = 1
   istride(3) = 1
   icount(1) = 1
   icount(2) = 1
   icount(3) = 1
   iblock(1) = dimsf_loc(1)
   iblock(2) = dimsf_loc(2)
   iblock(3) = dimsf_loc(3)

   if(myid .eq. 0) then
      dimsf_loc(1)=nx_loc+2*bufx
      dimsf_loc(2)=ny_loc+2*bufy
      dimsf_loc(3)=nz_loc+bufz
      istart(1) = 0
      istart(2) = 0      
      istart(3) = 0
      iblock(1) = dimsf_loc(1)
      iblock(2) = dimsf_loc(2)
      iblock(3) = dimsf_loc(3)
      allocate(q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc))
      q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc)=q_loc(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc)
   endif

   if(myid .eq. nprocs-1) then
      dimsf_loc(1)=nx_loc+2*bufx
      dimsf_loc(2)=ny_loc+2*bufy
      dimsf_loc(3)=nz_loc+bufz
      iblock(1) = dimsf_loc(1)
      iblock(2) = dimsf_loc(2)
      iblock(3) = dimsf_loc(3)
      allocate(q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1:nz_loc+bufz))
      q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1:nz_loc+bufz)=q_loc(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1:nz_loc+bufz)
   endif

   if(myid .ne. 0 .and. myid .ne. nprocs-1) then
      allocate(q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1:nz_loc))
      q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1:nz_loc)=q_loc(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1:nz_loc)
   endif

else
   dimsf_loc(1) = dimsf(1)
   dimsf_loc(2) = dimsf(2)
   dimsf_loc(3) = dimsf(3)
   istart(1)=0
   istart(2)=0
   istart(3)=0
   istride(1)=1
   istride(2)=1
   istride(3)=1
   icount(1)=1
   icount(2)=1
   icount(3)=1
   iblock(1)=dimsf_loc(1)
   iblock(2)=dimsf_loc(2)
   iblock(3)=dimsf_loc(3)
   allocate(q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc+bufz))
   q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc+bufz)=q_loc(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc+bufz)
endif

!! create file/memory space
call h5screate_simple_f(3,dimsf,filespace,error)
call h5screate_simple_f(3,dimsf_loc,mspace_id,error)
!! create dataset
call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,istart,icount,error,istride,iblock)

call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace, dset_id,error)


call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,q_temp,dimsf_loc,error, file_space_id=filespace,mem_space_id=mspace_id,xfer_prp=plist_id)
call h5dclose_f(dset_id,error)
call h5pclose_f(plist_id,error)

call h5sclose_f(filespace,error)
call h5sclose_f(mspace_id,error)
call h5fclose_f(file_id,error)
call h5close_f(error)

deallocate(q_temp)

end subroutine
