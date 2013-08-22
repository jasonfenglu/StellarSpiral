PROGRAM rotation
USE STELLARDISK_MODEL
USE STELLARDISK
IMPLICIT NONE
type(typspiral)                         ::spiral
CHARACTER(len=40),PARAMETER             ::datfname='data/RC_Lowe'
INTEGER,PARAMETER                       ::datlength=20
DOUBLE PRECISION                        ::RcDat(datlength,2)
DOUBLE PRECISION                        ::r,dr
DOUBLE PRECISION,ALLOCATABLE            ::dat(:,:)
INTEGER                                 ::i
INTEGER,PARAMETER                       ::N = 1000
DOUBLE PRECISION,PARAMETER              ::rf = 25.d0,ri = 0.d0


!Read Kendall's data into dat first 4 rows
open(10,file=datfname,ACTION='READ')
DO i = 1, datlength
        read(10,*)Rcdat(i,1:2)
ENDDO
close(10)

dr = (rf-ri)/dble(N)
ALLOCATE(dat(N,6))
CALL stdpara.readstd
CALL spiral.readw(2)
CALL spiral.init(500,10.d0,stdpara,2)
DO i = 1, N
        r = dr*dble(i)
        dat(i,1) = r
        dat(i,2) = r*Omega(r,spiral)
        dat(i,3) = r*bOmega(r,spiral)
        dat(i,4) = r*dOmega(r,spiral)
        dat(i,5) = r*hOmega(r,spiral)
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
CALL PGENV(0.,25.,0.,250.,0,0)

CALL PGPT(datlength,real(RcDat(:,1)),real(RcDat(:,2)),2)

DO i = 2, 5
        CALL PGSCI(i)
        CALL PGLINE(N,real(dat(:,1)),real(dat(:,i)))
ENDDO

CALL PGSCI(1)

CALL PGSCI(1)
CALL PGLAB('Radius','Rotation Velocity','')
CALL PGCLOS

ENDSUBROUTINE

ENDPROGRAM
