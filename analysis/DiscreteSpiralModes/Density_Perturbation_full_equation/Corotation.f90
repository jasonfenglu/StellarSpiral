PROGRAM corotation
USE STELLARDISK_MODEL
USE STELLARDISK
IMPLICIT NONE
type(typspiral)                         ::spiral
DOUBLE PRECISION,ALLOCATABLE            ::dat(:,:)
DOUBLE PRECISION,PARAMETER              ::rf = 12.d0,ri = 0.d0
DOUBLE PRECISION                        ::dr,k,O,r
INTEGER,PARAMETER                       ::N = 1000
INTEGER                                 ::i

dr = (rf-ri)/dble(N)
CALL stdpara.readstd
CALL spiral.readw(2)
CALL spiral.init(500,12.d0,stdpara,2)
CALL FindSpiral(spiral)
ALLOCATE(dat(N,7))
DO i = 1, N
        r = dr*dble(i)
        k = kappa(r,spiral)
        o = Omega(r,spiral)
        dat(i,1) =  r
        dat(i,2) =  o
        dat(i,3) =  o - k/2.d0
        dat(i,4) =  o + k/2.d0
        dat(i,5) =  o - k/4.d0
        dat(i,6) =  o + k/4.d0
        dat(i,7) =  real(spiral.w)/2.d0
ENDDO

CALL plot(dat)
CONTAINS

SUBROUTINE plot(dat)
DOUBLE PRECISION,INTENT(IN)             ::dat(:,:)
CHARACTER(len=40)                       ::text
INTEGER                                 ::PGBEG
INTEGER                                 ::i,j
IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
CALL PGSVP(0.0,0.95,0.0,0.95)
CALL PGENV(0.,12.,0.,120.,0,0)
DO i = 2, 7
        CALL PGSCI(i)
        CALL PGLINE(N,real(dat(:,1)),real(dat(:,i)))
ENDDO

CALL PGSCI(1)
write(text,"('Pattern Speed:',G12.4)")real(spiral.w)/2.d0
CALL PGPTXT(6.,100.,0.,0.,text)
write(text,"('Corotation:',G12.4)")CO()
CALL PGPTXT(6.,90.,0.,0.,text)
write(text,"('4 to 1:',G12.4)")four21()
CALL PGPTXT(6.,80.,0.,0.,text)

CALL PGSCI(1)
CALL PGLAB('Radius','Rotation Velocity','')
CALL PGCLOS

ENDSUBROUTINE

FUNCTION CO()
DOUBLE PRECISION                        ::CO
DOUBLE PRECISION                        ::B,C,R,RE,AE
INTEGER                                 ::IFLAG

B = 5.d0
C = 14.d0
R = 7.d0
RE = 1d-7
AE = 1d-7

CALL DFZERO(fco,B,C,R,RE,AE,IFLAG)
CO = c
ENDFUNCTION

FUNCTION four21()
DOUBLE PRECISION                        ::four21
DOUBLE PRECISION                        ::B,C,R,RE,AE
INTEGER                                 ::IFLAG

B = 5.d0
C = 14.d0
R = 7.d0
RE = 1d-7
AE = 1d-7

CALL DFZERO(ffour21,B,C,R,RE,AE,IFLAG)
four21 = c
ENDFUNCTION

FUNCTION fco(r)
DOUBLE PRECISION                ::fco,r
        fco = Omega(r,spiral)-real(spiral.w)/2.d0
ENDFUNCTION

FUNCTION ffour21(r)
DOUBLE PRECISION                ::ffour21,r
        ffour21 = Omega(r,spiral) + kappa(r,spiral)/4.d0 -real(spiral.w)/2.d0
ENDFUNCTION
ENDPROGRAM
