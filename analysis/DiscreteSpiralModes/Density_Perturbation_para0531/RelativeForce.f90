PROGRAM kendall
USE STELLARDISK_MODEL
USE STELLARDISK
USE io
IMPLICIT NONE
CHARACTER(len=40),PARAMETER             ::forceh5nm='force.h5'
type(typspiral)                         ::spiral
DOUBLE PRECISION,ALLOCATABLE            ::f(:,:,:)
DOUBLE PRECISION,ALLOCATABLE            ::x(:),y(:)
DOUBLE PRECISION,ALLOCATABLE            ::lf(:,:)
DOUBLE PRECISION                        ::r
INTEGER                                 ::fsize(2)
INTEGER                                 ::i,j

CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)

!read force map from force.h5
CALL readforce

ALLOCATE(lf(fsize(1)*fsize(2),2))
DO I = 1, fsize(1)
DO J = 1, fsize(2)
        r = sqrt(x(I)**2+y(J)**2)
        lf(I+(J-1)*fsize(1),1) = r
        lf(I+(J-1)*fsize(1),2) = sqrt(f(i,j,2)**2 + f(i,j,2)**2) &
                               /Omega(r,spiral)**2/r
ENDDO
ENDDO
CALL plot

CONTAINS

SUBROUTINE readforce
CHARACTER(len=40),PARAMETER             ::forceh5nm='force.h5'
DOUBLE PRECISION,ALLOCATABLE            ::tmp(:)
INTEGER,PARAMETER                       ::buf = 2
INTEGER                                 ::ndim(2)
ndim = h5size(forceh5nm,"density",2)

ALLOCATE(f(ndim(1),ndim(2),4))
CALL h5read(f(:,:,1),ndim(1),ndim(2),forceh5nm,"density")
CALL h5read(f(:,:,2),ndim(1),ndim(2),forceh5nm,"stellarfx")
CALL h5read(f(:,:,3),ndim(1),ndim(2),forceh5nm,"stellarfy")

ALLOCATE(tmp(ndim(1)+2*buf))
ALLOCATE(x(ndim(1)))
ALLOCATE(y(ndim(1)))
CALL h5read(tmp,size(tmp),forceh5nm,"x")
x = tmp(1+buf:ndim(1)+buf)
CALL h5read(tmp,size(tmp),forceh5nm,"y")
y = tmp(1+buf:ndim(1)+buf)

fsize = ndim

ENDSUBROUTINE

SUBROUTINE plot
INTEGER                                 ::PGBEG
INTEGER                                 ::i,j
IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
CALL output
!IF (PGBEG(0,'RelativeForce.ps/ps',1,1) .NE. 1) STOP
!CALL output

ENDSUBROUTINE

SUBROUTINE output
REAL                                    ::TR(6)         !plot geometry
CALL PGSVP(0.0,0.95,0.0,0.95)
CALL PGENV(0.,real(maxval(x))*sqrt(2.),0.,0.5,0,0)
!CALL PGENV(7.,8.,0.,0.1,0,0)

CALL PGPT(fsize(1)*fsize(2),real(lf(:,1)),real(lf(:,2)),1)

CALL PGSCI(1)
CALL PGCLOS

ENDSUBROUTINE

ENDPROGRAM

