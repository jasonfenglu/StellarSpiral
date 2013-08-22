MODULE io
USE HDF5
IMPLICIT NONE
TYPE::tyfiles
        CHARACTER(LEN=20)       ::filename
        TYPE(tyfiles),POINTER   ::ptnext
ENDTYPE
TYPE(tyfiles),POINTER,SAVE      ::ptfiles
TYPE(tyfiles),POINTER           ::ptfile

INTERFACE h5write
        MODULE PROCEDURE w2d
        MODULE PROCEDURE w1d
ENDINTERFACE

INTERFACE h5size
        MODULE PROCEDURE h5size1d
        MODULE PROCEDURE h5size2d
ENDINTERFACE

INTERFACE h5read
        MODULE PROCEDURE r1d
        MODULE PROCEDURE r2d
        MODULE PROCEDURE r1da
        MODULE PROCEDURE r2da
ENDINTERFACE

CONTAINS

FUNCTION h5size1d(filename,dsetname,ndim)
IMPLICIT NONE
CHARACTER(LEN=*)                :: filename      ! File name
CHARACTER(LEN=*)                :: dsetname      ! Dataset name
INTEGER(HID_T)                  :: file_id       ! File identifier
INTEGER(HID_T)                  :: dset_id       ! Dataset identifier
INTEGER(HID_T)                  :: dspace_id     ! Dataspace identifier
INTEGER(HSIZE_T)                :: dims(ndim)       ! Dataset dimensions
INTEGER(HSIZE_T)                :: maxdims(ndim)    ! Max dataset dimensions
INTEGER                         :: hdferr        ! Error flag
INTEGER                         :: ndim
INTEGER                         :: h5size1d(ndim)

CALL h5open_f(hdferr)
CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)
CALL h5dopen_f(file_id, dsetname, dset_id, hdferr)
CALL h5dget_space_f(dset_id, dspace_id, hdferr) 
CALL h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr) 
CALL h5dclose_f(dset_id, hdferr)
CALL h5sclose_f(dspace_id, hdferr)
CALL h5fclose_f(file_id, hdferr)
CALL h5close_f(hdferr)

h5size1d = int(dims)

ENDFUNCTION

FUNCTION h5size2d(filename,dsetname)
IMPLICIT NONE
CHARACTER(LEN=*)                :: filename      ! File name
CHARACTER(LEN=*)                :: dsetname      ! Dataset name
INTEGER(HID_T)                  :: file_id       ! File identifier
INTEGER(HID_T)                  :: dset_id       ! Dataset identifier
INTEGER(HID_T)                  :: dspace_id     ! Dataspace identifier
INTEGER(HSIZE_T)                :: dims(2)       ! Dataset dimensions
INTEGER(HSIZE_T)                :: maxdims(2)    ! Max dataset dimensions
INTEGER                         :: hdferr        ! Error flag
INTEGER                         :: h5size2d(2)

CALL h5open_f(hdferr)
CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)
CALL h5dopen_f(file_id, dsetname, dset_id, hdferr)
CALL h5dget_space_f(dset_id, dspace_id, hdferr) 
CALL h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr) 
CALL h5dclose_f(dset_id, hdferr)
CALL h5sclose_f(dspace_id, hdferr)
CALL h5fclose_f(file_id, hdferr)
CALL h5close_f(hdferr)

h5size2d = int(dims)

ENDFUNCTION

SUBROUTINE r1d(dat,M,filename,dsetname)
CHARACTER(LEN=*)                :: filename      ! File name
CHARACTER(LEN=*)                :: dsetname      ! Dataset name
DOUBLE PRECISION,INTENT(OUT)    :: dat(:)        ! input data
INTEGER(HID_T)                  :: file_id       ! File identifier
INTEGER(HID_T)                  :: dset_id       ! Dataset identifier
INTEGER(HSIZE_T)                :: dims(1)       ! Dataset dimensions
INTEGER                         :: rank = 1      ! Dataset rank
INTEGER                         :: error         ! Error flag
INTEGER                         :: M             ! Dimension of data

dims = M

CALL h5open_f(error)
CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
CALL h5dopen_f(file_id, dsetname, dset_id, error)
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat, dims, error)
CALL h5dclose_f(dset_id, error)
CALL h5fclose_f(file_id, error)
CALL h5close_f(error)

ENDSUBROUTINE

SUBROUTINE r1da(dat,filename,dsetname)
CHARACTER(LEN=*)                :: filename      ! File name
CHARACTER(LEN=*)                :: dsetname      ! Dataset name
DOUBLE PRECISION,INTENT(OUT),ALLOCATABLE:: dat(:)! input data
INTEGER(HID_T)                  :: file_id       ! File identifier
INTEGER(HID_T)                  :: dset_id       ! Dataset identifier
INTEGER(HSIZE_T)                :: dims(1)       ! Dataset dimensions
INTEGER                         :: rank = 1      ! Dataset rank
INTEGER                         :: error         ! Error flag

dims = h5size(filename,dsetname,1)
ALLOCATE(dat(dims(1)))

CALL h5open_f(error)
CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
CALL h5dopen_f(file_id, dsetname, dset_id, error)
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat, dims, error)
CALL h5dclose_f(dset_id, error)
CALL h5fclose_f(file_id, error)
CALL h5close_f(error)

ENDSUBROUTINE

SUBROUTINE r2d(dat,M,N,filename,dsetname)
CHARACTER(LEN=*)                :: filename      ! File name
CHARACTER(LEN=*)                :: dsetname      ! Dataset name
DOUBLE PRECISION,INTENT(OUT)    :: dat(:,:)      ! input data
INTEGER(HID_T)                  :: file_id       ! File identifier
INTEGER(HID_T)                  :: dset_id       ! Dataset identifier
INTEGER(HSIZE_T)                :: dims(2)       ! Dataset dimensions
INTEGER                         :: rank = 2      ! Dataset rank
INTEGER                         :: error         ! Error flag
INTEGER                         :: M,N           ! Dimension of data

dims = (/M,N/)

CALL h5open_f(error)
CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
CALL h5dopen_f(file_id, dsetname, dset_id, error)
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat, dims, error)
CALL h5dclose_f(dset_id, error)
CALL h5fclose_f(file_id, error)
CALL h5close_f(error)

ENDSUBROUTINE

SUBROUTINE r2da(dat,filename,dsetname)
CHARACTER(LEN=*)                :: filename      ! File name
CHARACTER(LEN=*)                :: dsetname      ! Dataset name
DOUBLE PRECISION,INTENT(OUT),ALLOCATABLE:: dat(:,:)! input data
INTEGER(HID_T)                  :: file_id       ! File identifier
INTEGER(HID_T)                  :: dset_id       ! Dataset identifier
INTEGER(HSIZE_T)                :: dims(2)       ! Dataset dimensions
INTEGER                         :: rank = 2      ! Dataset rank
INTEGER                         :: error         ! Error flag

dims = h5size(filename,dsetname)

if(.not.ALLOCATED(dat))ALLOCATE(dat(dims(1),dims(2)))

CALL h5open_f(error)
CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
CALL h5dopen_f(file_id, dsetname, dset_id, error)
CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat, dims, error)
CALL h5dclose_f(dset_id, error)
CALL h5fclose_f(file_id, error)
CALL h5close_f(error)

ENDSUBROUTINE

SUBROUTINE w2d(dat,M,N,filename,dsetname)
CHARACTER(LEN=*)                :: filename      ! File name
CHARACTER(LEN=*)                :: dsetname      ! Dataset name
DOUBLE PRECISION,INTENT(IN)     :: dat(:,:)
INTEGER(HID_T)                  :: file_id       ! File identifier
INTEGER(HID_T)                  :: dset_id       ! Dataset identifier
INTEGER(HID_T)                  :: dspace_id     ! Dataspace identifier
INTEGER(HSIZE_T)                :: dims(2)     ! Dataset dimensions
INTEGER                         :: rank = 2      ! Dataset rank
INTEGER                         :: error         ! Error flag
INTEGER                         :: M,N           ! Dimension of data

dims = (/M,N/)


CALL h5open_f(error)
IF(checkfileext(filename))THEN     
        CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
ELSE
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
ENDIF
!prepare dataspace
CALL h5screate_simple_f(rank, dims, dspace_id, error)
!create dataset
CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
     dset_id, error)
!! write the dataset
CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat, dims, error)
!closing
CALL h5dclose_f(dset_id, error)
CALL h5sclose_f(dspace_id, error)
CALL h5fclose_f(file_id, error)
CALL h5close_f(error)

ENDSUBROUTINE

SUBROUTINE w1d(dat,M,filename,dsetname)
CHARACTER(LEN=*)                :: filename      ! File name
CHARACTER(LEN=*)                :: dsetname      ! Dataset name
DOUBLE PRECISION,INTENT(IN)     :: dat(:)        ! input data
INTEGER(HID_T)                  :: file_id       ! File identifier
INTEGER(HID_T)                  :: dset_id       ! Dataset identifier
INTEGER(HID_T)                  :: dspace_id     ! Dataspace identifier
INTEGER(HSIZE_T)                :: dims(1)     ! Dataset dimensions
INTEGER                         :: rank = 1      ! Dataset rank
INTEGER                         :: error         ! Error flag
INTEGER                         :: M             ! Dimension of data

dims = M

CALL h5open_f(error)
IF(checkfileext(filename))THEN     
        CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
ELSE
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
ENDIF
!prepare dataspace
CALL h5screate_simple_f(rank, dims, dspace_id, error)
!create dataset
CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
     dset_id, error)
!! write the dataset
CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat, dims, error)
!closing
CALL h5dclose_f(dset_id, error)
CALL h5sclose_f(dspace_id, error)
CALL h5fclose_f(file_id, error)
CALL h5close_f(error)

ENDSUBROUTINE

FUNCTION checkfileext(fname)
CHARACTER(LEN=*),INTENT(IN)     ::fname
LOGICAL                         ::checkfileext
INTEGER                         ::err
IF(.not.ASSOCIATED(ptfiles))then
        !!if file list is not exist, create one and return false
        !print *,'file list not exist, create one.'
        ALLOCATE(ptfiles)
        checkfileext = .false.
        ptfiles.filename = fname
        ptfiles.ptnext => null()
        RETURN
ELSE
        !!file list exist, check if filename in the list
        ptfile => ptfiles
        IF(ptfile.filename==fname)then
                !print *,'filename in the list','first'
                checkfileext = .true.
                RETURN
        ENDIF
        DO WHILE(ASSOCIATED(ptfile.ptnext))
        !print *,'loop search'
        if(trim(ptfile.filename)==trim(fname))then
                !!fname is in the list
                !print *,'filename in the list'
                checkfileext = .true.
                RETURN
        ELSE
        ptfile => ptfile.ptnext
        ENDIF
        ENDDO
        !!file not in the list, create filename in the next obj
        !print *,'filename not in the list, create one',trim(ptfile.filename)==trim(fname)
        ALLOCATE(ptfile.ptnext,stat=err)
        ptfile => ptfile.ptnext
        ptfile.filename = fname
        ptfile.ptnext => null()
        checkfileext = .false.
        RETURN
ENDIF
ENDFUNCTION

SUBROUTINE printpt2
print *,'!!!!'
!ptfile => ptfiles
!print *,'@@@',ptfile.filename
!ptfile => ptfile.ptnext
!print *,'@@@',ptfile.filename
print *,'1',ptfiles.filename
print *,'2',ptfiles.ptnext.filename
print *,'!!!!'
ENDSUBROUTINE

SUBROUTINE printpt
INTEGER i 
print *,'!!!!','print stored filenaes'
!!check if file list exist
IF(.not.ASSOCIATED(ptfiles))then
        print *,'file list not exist'
        return
ENDIF
!!file list exist, print all

ptfile => ptfiles
print *,ptfile.filename
DO WHILE(ASSOCIATED(ptfile.ptnext))
        ptfile => ptfile.ptnext
        print *,ptfile.filename
ENDDO
print *,'!!!!'
ENDSUBROUTINE

ENDMODULE
