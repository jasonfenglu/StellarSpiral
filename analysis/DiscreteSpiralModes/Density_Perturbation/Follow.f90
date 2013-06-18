MODULE FollowMO
USE math,only:typintplt2
IMPLICIT NONE
TYPE    tyH5file
        CHARACTER(len=255)              ::filename
        DOUBLE PRECISION,ALLOCATABLE    ::density(:,:)
        DOUBLE PRECISION,ALLOCATABLE    ::x(:)
        DOUBLE PRECISION,ALLOCATABLE    ::y(:)
        INTEGER                         ::n
        CONTAINS
        PROCEDURE,PASS::init            =>new_tyH5file
        PROCEDURE,PASS::free            =>free_tyH5file
ENDTYPE
TYPE(tyH5file)                          ::gas

TYPE    tyFollower
        DOUBLE PRECISION                ::x,y,r,th,amp
        CONTAINS
        PROCEDURE,PASS::init            =>new_Follower
ENDTYPE

!!2d intepolation 
type(typintplt2)                        ::intplt2
CONTAINS

SUBROUTINE new_Follower(this,x,y)
IMPLICIT NONE
CLASS(tyFollower)                       ::this
DOUBLE PRECISION                        ::x,y

this.x = x
this.y = y
this.r = sqrt(this.x**2+this.y**2)
this.th= atan2(this.y,this.x)
this.amp=intplt2.find(this.x,this.y)

ENDSUBROUTINE new_Follower

SUBROUTINE new_tyH5file(this,filename)
USE io
IMPLICIT NONE
CLASS(tyH5file)                         ::this
CHARACTER(len=255)                      ::filename
this.filename = filename
CALL h5read(this.density,this.filename,'density')
CALL h5read(this.x,this.filename,'x')
CALL h5read(this.y,this.filename,'y')
this.n = size(this.x)
ENDSUBROUTINE new_tyH5file

SUBROUTINE free_tyH5file(this)
USE io
IMPLICIT NONE
CLASS(tyH5file)                         ::this
DEALLOCATE(this.density)
DEALLOCATE(this.x)
DEALLOCATE(this.y)
ENDSUBROUTINE free_tyH5file

ENDMODULE

MODULE argument
USE FollowMO
IMPLICIT  NONE
        CHARACTER(len=255)      ::filename
        DOUBLE PRECISION        ::zmax      = 30.d0
        INTEGER                 ::xn
        LOGICAL                 ::xn_set    = .false.
        LOGICAL                 ::zmax_set  = .false.
        LOGICAL                 ::givefnm   = .false.

CONTAINS

SUBROUTINE readarg
IMPLICIT NONE
CHARACTER(len=32)               ::arg
INTEGER                         ::i
INTEGER                         ::ioerr
!read in rotating parameters and gas filename by stdin
if(iargc().ne.0)then
        DO i = 1, iargc()
                CALL getarg(i,arg)
                SELECT CASE(arg)
                CASE('-i')
                        CALL getarg(i+1,arg)
                        READ(arg,*)filename
                        givefnm = .true.
                CASE('--help','-h')
                        CALL help
                CASE('--zmax','-z')
                        CALL getarg(i+1,arg)
                        READ(arg,*)zmax
                        zmax_set = .true.
                CASE('-x')
                        CALL getarg(i+1,arg)
                        READ(arg,*)xn
                        xn_set = .true.
                CASE('-v')
                        write(6,'(a)')'Compiled at: '//__DATE__//' '//__TIME__
                        STOP
                ENDSELECT
        ENDDO
ENDIF
ENDSUBROUTINE

SUBROUTINE help
IMPLICIT NONE
        write(6,'(a)')'Follow spiral.                       '
        write(6,'(a)')'usage:   Follow.exe [options]               '
        write(6,'(a)')'options:                         '
        write(6,'(a)')'         -z, --zmax [scale of z] Specified the scale of z'
        write(6,'(a)')'         -h, --help              Show this help page     '
        write(6,'(a)')'         -v,                     Version information.'
        write(6,'(a)')'         -x,                     Number of XWindows. '
        STOP
ENDSUBROUTINE

ENDMODULE

PROGRAM FOLLOW
USE argument
IMPLICIT NONE
TYPE(tyFollower)                        ::tmp,tmp2
INTEGER,PARAMETER                       ::N = 150
INTEGER,PARAMETER                       ::M = 2
TYPE(tyFollower)                        ::spiral(N,M)
TYPE(tyFollower)                        ::FirstFollower
INTEGER                                 ::i,j

!!reading from arguments
CALL readarg

!!reading data
CALL gas.init(filename)

!!init intepolating obj
CALL intplt2.init(gas.density,gas.x,gas.y)

!!first plot
CALL InitPlot

!!plot circle
!CALL PGSAVE
!CALL PGSCI(0)
!CALL PGSFS(2)
!CALL PGCIRC(0.,0.,3.)
!CALL PGCIRC(0.,0.,4.)
!CALL PGCIRC(0.,0.,5.)
!CALL PGUNSA

!!read by click
DO j = 1, M
CALL FirstPoint(FirstFollower)
CALL NextPoint(FirstFollower,spiral(1,j))
DO i = 2, N
        CALL NextPoint(spiral(i-1,j),spiral(i,j))
ENDDO
CALL PGSAVE
CALL PGSCI(j+1)
CALL PGPT(N,real(spiral(:,j).x),real(spiral(:,j).y),2)
CALL PGUNSA
ENDDO
CALL PlotSpiral(N,M,spiral)

!!ending program
CALL gas.free
CALL PGCLOS
STOP
ENDPROGRAM

SUBROUTINE InitPlot
USE FollowMO
USE plotting
USE argument
IMPLICIT NONE
INTEGER                                 ::PGOPEN
DOUBLE PRECISION                        ::domain
IF(.not.zmax_set)zmax = maxval(gas.density)
write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'

IF(xn_set)THEN
        IF(PGOPEN(ACHAR(48+xn)//'/xs').NE.1)STOP
ELSE
        IF (PGOPEN('/xs') .NE. 1) STOP
ENDIF

CALL PGSUBP(2,1)
domain = maxval(gas.x)
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
CALL meshplot(gas.density,gas.n,domain,zmax,zmin_in=0.d0)
ENDSUBROUTINE

SUBROUTINE PlotSpiral(N,M,Followers)
USE FollowMO
USE plotting
USE argument
IMPLICIT NONE
INTEGER                                         ::N,M
TYPE(tyFollower)                                ::Followers(N,M)
DOUBLE PRECISION                                ::thi,thf,ampi,ampf
DOUBLE PRECISION                                ::ri,rf
INTEGER                                         ::i
thi = minval(Followers(:,:).th)*0.9
thf = maxval(Followers(:,:).th)*1.1
ampi = minval(Followers(:,:).amp)*0.9
ampf = maxval(Followers(:,:).amp)*1.1
ri  = minval(Followers(:,:).r)*0.9
rf  = maxval(Followers(:,:).r)*1.1

ampi = 20.
IF(zmax_set)ampf = zmax

CALL PGPAGE
!!plot th to amp
CALL PGSWIN(real(thf),real(thi),real(ampi),real(ampf))
CALL PGSVP(0.1,0.9,0.6,0.9)
CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
CALL PGSAVE
DO i = 1, M
        CALL PGSCI(i+1)
        CALL PGLINE(N,real(Followers(:,i).th),real(Followers(:,i).amp))
ENDDO
CALL PGUNSA
CALL PGLAB('\gh','Amp','\gh to Amp')

!!plot th to amp
CALL PGSWIN(real(ri),real(rf),real(ampi),real(ampf))
CALL PGSVP(0.1,0.9,0.1,0.4)
CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
CALL PGSAVE
DO i = 1, M
        CALL PGSCI(i+1)
        CALL PGLINE(N,real(Followers(:,i).r),real(Followers(:,i).amp))
ENDDO
CALL PGUNSA
CALL PGLAB('r','Amp','r to Amp')
ENDSUBROUTINE

SUBROUTINE FirstPoint(Follower)
USE FollowMO
IMPLICIT NONE
TYPE(tyFollower)                ::Follower
CHARACTER(1)                    ::CH
DOUBLE PRECISION                ::PGCURS
REAL                            ::rx,ry
INTEGER                         ::ierr

ierr = PGCURS(rx,ry,CH)
CALL Follower.init(dble(rx),dble(ry))

ENDSUBROUTINE

SUBROUTINE NextPoint(Follower,NxFollower)
USE FollowMO
USE math
IMPLICIT NONE
TYPE(tyFollower),INTENT(IN)     ::Follower
TYPE(tyFollower),INTENT(OUT)    ::NxFollower
DOUBLE PRECISION,PARAMETER      ::r=0.1d0
DOUBLE PRECISION                ::th,dth,amp
DOUBLE PRECISION                ::pi
INTEGER,PARAMETER               ::N = 600
DOUBLE PRECISION                ::dat(N,2)
INTEGER                         ::i


pi = atan(1.d0)*4.d0
dth = -pi/dble(N) 
DO i = 1, N
        th = dth*dble(i)+ atan2(Follower.y,Follower.x)
        amp= intplt2.find(Follower.x+r*cos(th),Follower.y+r*sin(th))
        dat(i,1) = th
        dat(i,2) = amp
ENDDO
th = dat(maxloc(dat(:,2),1),1)
CALL NxFollower.init(Follower.x+r*cos(th),Follower.y+r*sin(th))

ENDSUBROUTINE

