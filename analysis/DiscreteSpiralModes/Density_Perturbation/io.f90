MODULE io
USE HDF5
IMPLICIT NONE
TYPE::tyfiles
        CHARACTER(LEN=20)       ::filename
        TYPE(tyfiles),POINTER   ::ptnext
ENDTYPE
TYPE(tyfiles),POINTER,SAVE      ::ptfiles
TYPE(tyfiles),POINTER           ::ptfile

INTERFACE h5io
        MODULE PROCEDURE io2d
        MODULE PROCEDURE io1d
ENDINTERFACE
CONTAINS

SUBROUTINE io2d(dat,M,N,filename,dsetname)
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

SUBROUTINE io1d(dat,M,filename,dsetname)
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
