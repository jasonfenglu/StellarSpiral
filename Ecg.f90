MODULE ECGMO
        CHARACTER(len=225)              ::fname
        DOUBLE PRECISION,ALLOCATABLE    ::density(:,:)
        DOUBLE PRECISION,ALLOCATABLE    ::x(:),y(:)
        DOUBLE PRECISION,SAVE           ::zmax
        LOGICAL                         ::zmax_set = .false.
ENDMODULE 

PROGRAM ECG
USE io
USE ECGMO
USE math
USE STELLARDISK_MODEL
USE STELLARDISK,only:FindSpiral,sigma1
IMPLICIT NONE
CHARACTER(len=32)               ::arg
type(typintplt2)                ::intplt2
type(typspiral)                 ::spiral
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
                CASE('--zmax','-z')
                        CALL getarg(i+1,arg)
                        READ(arg,*)zmax
                        zmax_set = .true.
                CASE('-v')
                        write(6,'(a)')'Compiled at: '//__DATE__//' '//__TIME__
                        STOP
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

!!init stellar spiral
CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)

!!set coordinate
ri = 0.d0
rf = spiral.co
dr = (rf-ri)/dble(N-1)
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
                dat(j,i) =  intplt2.find(r*cos(th),r*sin(th))
!               dat(j,i) =  sigma1(r,th,spiral)
        ENDDO
ENDDO

CALL plot(density,size(x),maxval(x))
CALL PGCLOS

STOP
CONTAINS
SUBROUTINE plot(F,n,domain)
USE plotting
USE ECGMO
IMPLICIT NONE
CHARACTER(len=32)                       ::text
DOUBLE PRECISION,INTENT(IN)             ::F(:,:)        !plotting data
DOUBLE PRECISION,INTENT(IN)             ::domain        !plot range
INTEGER,INTENT(IN)                      ::n             !dimentsion
INTEGER                                 ::PGBEG
INTEGER                                 ::i

IF(.not.zmax_set)zmax = maxval(F)
write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'

IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP

CALL PGSUBP(2,1)
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
CALL meshplot(F,n,domain,zmax,zmin_in=0.d0)
CALL PGSFS(2)
CALL PGSCI(6)
DO i = 1, 5
        CALL PGCIRC(0.,0.,real(rs(i)))
ENDDO
CALL PGSCI(1)

!!move to right panel
CALL PGPAGE
CALL PGSWIN(0.0,real(3.d0*pi),0.,80.)
!CALL PGSWIN(0.0,real(3.d0*pi),-800.,800.)
DO i = 2, 5
        CALL PGSAVE
        CALL PGSVP(0.1,0.9,0.1+real(i-2)*0.2,0.25+real(i-2)*0.2)
        print *,0.25+real(i-2)*0.2
        CALL PGBOX('ABCTSN',0.0,0,'ABCTSN',0.0,0)
        CALL PGSCI(i)
        CALL PGLINE(M,real(dat(:,1)),real(dat(:,i)))
        write(text,'(G11.4)')rs(i)
        CALL PGTEXT(1.,60.,text)
        CALL PGUNSA
ENDDO

ENDSUBROUTINE
ENDPROGRAM ECG

SUBROUTINE help
write(6,'(a)')'Show amplitude of some radius.                   '
write(6,'(a)')'usage:   Ecg.exe -i filename.h5 [-z zmax]        '
write(6,'(a)')'         Ecg.exe -v                              '
STOP
ENDSUBROUTINE help
