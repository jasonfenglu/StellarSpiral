PROGRAM kendall
USE STELLARDISK_MODEL
USE STELLARDISK
IMPLICIT NONE
CHARACTER(len=32)                       ::arg
CHARACTER(len=40),PARAMETER             ::datfname='data/ReAmpKendall.csv'
type(typspiral)                         ::spiral
INTEGER,PARAMETER                       ::datlength=207
DOUBLE PRECISION                        ::dat(datlength,7)
DOUBLE PRECISION                        ::Amp
DOUBLE PRECISION                        ::r,dr = 0.1d0
DOUBLE PRECISION                        ::bias = 1.d0
INTEGER                                 ::i

if(iargc().ne.0)then
        DO i = 1, iargc()
                CALL getarg(i,arg)
                SELECT CASE(arg)
                CASE('--help','-h')
                        write(6,'(a)')'Calculate prefactor of density contrast. '
                        write(6,'(a)')'usage:   Kendall.exe [option]            '
                        write(6,'(a)')'options:                         '
                        write(6,'(a)')'         -f, --force             Increase force multiple.'
                        STOP
                CASE('--force','-f')
                        CALL getarg(i+1,arg)
                        READ(arg,*)bias
                        print *,'using bias:',bias
                ENDSELECT
        ENDDO
ENDIF

!Read Kendall's data into dat first 4 rows
open(10,file=datfname,ACTION='READ')
read(10,*)
DO i = 1, datlength
        read(10,*)dat(i,1:4)
ENDDO
close(10)

!output to check
!DO i = 1, 207
!        write(6,'(4(G12.4,3X))')kdat(i,:)
!ENDDO

CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)

CALL findAmp(amp)
!amp = 1.d0
!!Fill in data
print *,'# Amp:',amp
DO i = 1, datlength
        r = dat(i,1) 
        dat(i,5) = amp*sigma1r(r,spiral)/(sigma0(r,spiral)+BulgeSurfaceDensity(r,spiral))
ENDDO

dat(:,6) = bias*dat(:,2)


CALL plot(dat)
CONTAINS

SUBROUTINE plot(dat)
DOUBLE PRECISION,INTENT(IN)             ::dat(:,:)
INTEGER                                 ::PGBEG
INTEGER                                 ::i,j
IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
CALL output
IF (PGBEG(0,'RelativeAmp.ps/ps',1,1) .NE. 1) STOP
CALL output

ENDSUBROUTINE

SUBROUTINE findAmp(amp)
USE STELLARDISK_MODEL
USE STELLARDISK
IMPLICIT NONE
!type(typspiral)                         ::spiral
!INTEGER                                 ::datlength
!DOUBLE PRECISION                        ::dat(:,:)
DOUBLE PRECISION                        ::amp
DOUBLE PRECISION                        ::B,C,R,RE,AE
INTEGER                                 ::IFLAG

B = 1.d3
C = 1.d8
R = 1.d6
RE = 1d-7
AE = 1d-7

CALL DFZERO(amperror,B,C,R,RE,AE,IFLAG)
amp = C

ENDSUBROUTINE

FUNCTION amperror(amp)
USE STELLARDISK_MODEL
USE STELLARDISK
IMPLICIT NONE
DOUBLE PRECISION                        ::amperror,amp
DOUBLE PRECISION                        ::r,tmp
INTEGER                                 ::i

amperror = 0.d0
DO i = 1, datlength
        r = dat(i,1)
        if(r.gt.spiral.fortoone)exit
        tmp = sigma1r(r,spiral)/(sigma0(r,spiral)+BulgeSurfaceDensity(r,spiral))
        amperror = amperror + (amp*tmp-bias*dat(i,2))**2
ENDDO
!print *,amperror,amp

ENDFUNCTION

SUBROUTINE output
CALL PGSVP(0.0,0.95,0.0,0.95)
CALL PGENV(3.,real(spiral.fortoone),0.,0.4,0,0)

DO i = 2, 6
        CALL PGSCI(i)
        !SELECT CASE(i)
        !CASE(2)
        !        CALL PGSLS(2)
        !CASE(3)
        !        CALL PGSLS(4)
        !CASE(4)
        !        CALL PGSLS(4)
        !CASE DEFAULT
        !        CALL PGSLS(1)
        !ENDSELECT
        CALL PGLINE(datlength,real(dat(:,1)),real(dat(:,i)))
        CALL PGSLS(1)
ENDDO

CALL PGSCI(1)
CALL PGLAB('Radius','Relative Ampltitude','')
CALL PGCLOS

ENDSUBROUTINE

ENDPROGRAM
