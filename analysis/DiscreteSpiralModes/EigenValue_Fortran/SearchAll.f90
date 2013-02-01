PROGRAM find_all
USE PLOTTING
USE STELLARDISK
IMPLICIT NONE
type searchgrid_type
        DOUBLE PRECISION::coord(4,4,2)
        DOUBLE PRECISION::error(4,4)
endtype
type(searchgrid_type)            ::searchgrid
DOUBLE PRECISION                  ::wri,wii
INTEGER                          ::l

DOUBLE PRECISION,ALLOCATABLE    ::a(:,:)

!!$OMP PARALLEL DEFAULT(NONE),PRIVATE(a)
!!$OMP DO
!DO l = 1,4
!        ALLOCATE(a(2,2))
!        print *,'jer',l
!        DEALLOCATE(a)
!ENDDO
!!$OMP END DO
!!$OMP END PARALLEL
!stop


wri = 63.d0
wii = -0.5d0
CALL findpspsd(wri,wii)
print *,wri,wii

STOP
END PROGRAM
