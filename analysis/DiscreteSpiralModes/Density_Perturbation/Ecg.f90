MODULE ECGMO
        CHARACTER(len=225)              ::fname
        DOUBLE PRECISION,ALLOCATABLE    ::density(:,:)
        DOUBLE PRECISION,ALLOCATABLE    ::x(:),y(:)
ENDMODULE 

PROGRAM ECG
USE io
USE ECGMO
USE math
IMPLICIT NONE
CHARACTER(len=32)               ::arg
type(typintplt2)                ::intplt2
DOUBLE PRECISION                ::r,th
INTEGER,PARAMETER               ::M = 1000
INTEGER,PARAMETER               ::N = 5         
DOUBLE PRECISION,PARAMETER      ::pi= 4.d0*atan(1.d0)
DOUBLE PRECISION                ::dat(M,N),dr,dth,ri,rf
DOUBLE PRECISION                ::rs(N)
INTEGER                         ::i,j,k,l
!Read arguments
IF(iargc().ne.0)THEN
        DO i = 1, iargc()
                CALL getarg(i,arg)
                SELECT CASE(arg)
                CASE('--help','-h')
                        CALL help
                CASE('--input','-i')
                        CALL getarg(i+1,arg)
                        READ(arg,*)fname
                ENDSELECT
        ENDDO
ELSE
        CALL help
ENDIF

!!reading h5 file
CALL h5read(density,fname,'density')
CALL h5read(x,fname,'x')
CALL h5read(y,fname,'y')

!!set interplating object
CALL intplt2.init(density,x,y)

!!set coordinate
ri = 1.d0
rf = 7.3d0
dr = (rf-ri)/dble(N)
dth = 3.d0*pi/dble(M)
DO j = 1, M
        th = dth * dble(j-1)
        dat(j,1) = th
ENDDO

DO i = 2, 5
        r = ri + dr*dble(i-1)
        rs(i) = r
        DO j = 1, M
                th = dth * dble(j-1)
                dat(j,i) =  intplt2.find(r*cos(th),r*sin(th)) + 50.d0*(i-2)
        ENDDO
ENDDO
CALL plot(density,size(x),maxval(x))
CALL PGCLOS

STOP
CONTAINS
SUBROUTINE plot(F,n,domain)
USE plotting
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)             ::F(:,:)        !plotting data
DOUBLE PRECISION,INTENT(IN)             ::domain        !plot range
INTEGER,INTENT(IN)                      ::n             !dimentsion
DOUBLE PRECISION                        ::zmax
INTEGER                                 ::PGBEG
INTEGER                                 ::i

zmax = maxval(F)*1.1d0
write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'

IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP

CALL PGSUBP(2,1)
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
CALL meshplot(F,n,domain,zmax,zmin_in=0.d0)
CALL PGSFS(2)
CALL PGSCI(0)
DO i = 1, 5
        CALL PGCIRC(0.,0.,real(rs(i)))
ENDDO
CALL PGSCI(1)

CALL PGENV(0.,real(maxval(dat(:,1))),0.,real(maxval(dat))*1.1,0,0)
DO i = 2, 5
        CALL PGSCI(i)
        CALL PGLINE(M,real(dat(:,1)),real(dat(:,i)))
ENDDO

ENDSUBROUTINE
ENDPROGRAM ECG

SUBROUTINE help
write(6,'(a)')'Print ECG like plot                              '
write(6,'(a)')'usage:   Ecg.exe -i filename.h5                  '
STOP
ENDSUBROUTINE help
