PROGRAM find_all
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

wri = 35.00
wii = -0.2
CALL findpspsd(wri,wii)
print *,wri,wii

STOP
END PROGRAM
