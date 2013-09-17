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
TYPE(tyH5file),TARGET                   ::gas

TYPE    tyFollower
        DOUBLE PRECISION                ::x,y,r,th,amp
        CONTAINS
        PROCEDURE,PASS::init            =>new_Follower
ENDTYPE

!!2d intepolation 
type(typintplt2),TARGET                 ::intpltgas
type(typintplt2),TARGET                 ::intpltstar
DOUBLE PRECISION,ALLOCATABLE            ::StellarDensity(:,:)
DOUBLE PRECISION                        ::pi = atan(1.d0)*4.d0
DOUBLE PRECISION                        ::GAC(58,2)
CONTAINS

SUBROUTINE new_Follower(this,x,y,intplt)
IMPLICIT NONE
CLASS(tyFollower)                       ::this
TYPE(typintplt2)                        ::intplt
DOUBLE PRECISION                        ::x,y

this.x = x
this.y = y
this.r = sqrt(this.x**2+this.y**2)
this.th= atan2(this.y,this.x)
this.amp=intplt.find(this.x,this.y)

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
        LOGICAL                 ::plot_circ = .false.
        LOGICAL                 ::output    = .false.

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
                CASE('-c')
                        plot_circ = .true.
                CASE('--output')
                        output = .true.
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
        write(6,'(a)')'Trace the amplitude of both stellar and gaseous spirals.'
        write(6,'(a)')'usage:   Follow.exe [options] -i Filename.h5'
        write(6,'(a)')'options:                         '
        write(6,'(a)')'         -z, --zmax [scale of z] Specified the scale of z'
        write(6,'(a)')'         -h, --help              Show this help page     '
        write(6,'(a)')'         -v,                     Version information.'
        write(6,'(a)')'         -x,                     Number of XWindows. '
        write(6,'(a)')'         -c,                     Plot circles at 3,4 and 5kpc. '
        write(6,'(a)')'         --output                Output Phase Shift. '
        STOP
ENDSUBROUTINE

ENDMODULE

PROGRAM Shift
USE argument
IMPLICIT NONE
TYPE(tyFollower)                        ::tmp,tmp2
INTEGER,PARAMETER                       ::N = 200
INTEGER,PARAMETER                       ::M = 2
TYPE(tyFollower)                        ::spiral(N,M)!!M=1 gas,M=2 star
INTEGER                                 ::i,j

!!reading from arguments
CALL readarg

!!reading data
CALL gas.init(filename)

!!init intepolating obj
CALL intpltgas.init(gas.density,gas.x,gas.y)

!!fill stellar plot
CALL InitStellar
CALL intpltstar.init(StellarDensity,gas.x,gas.y)

!!Plotting
CALL InitPlot
!!gas density plot
CALL InitGasPlot
!!star density plot
CALL InitStarPlot

!!Tracing spirals
write(6,'(a)')'Click initial tracing points of gas spiral:'
CALL TraceSpiral('gas',N,spiral(:,1))
write(6,'(a)')'Click initial tracing points of star spiral:'
CALL TraceSpiral('star',N,spiral(:,2))

!!Plot Results
CALL PlotSpiral(N,M,spiral)
CALL PGPANL(3,1)
CALL PlotPhase(N,M,spiral)
IF(output)CALL OutputPhase(N,M,spiral)

!!ending program
CALL gas.free
CALL FreeStellar
CALL PGCLOS
CALL intpltstar.free
CALL intpltgas.free
STOP
ENDPROGRAM

SUBROUTINE TraceSpiral(typ,N,follower)
USE FollowMO
USE math,only:typintplt2
IMPLICIT NONE
CHARACTER(*)                            ::typ
INTEGER                                 ::N
TYPE(tyFollower)                        ::follower(N)
TYPE(tyFollower)                        ::FirstFollower
type(typintplt2),POINTER                ::intplt=>null()
INTEGER                                 ::i,j

SELECT CASE(typ)
CASE('gas')
        intplt => intpltgas
        CALL PGPANL(1,1)
CASE('star')
        intplt => intpltstar
        CALL PGPANL(1,2)
ENDSELECT

!!hand input
CALL FirstPoint(FirstFollower,intplt)
!!correct hand input
CALL NextPoint(FirstFollower,follower,intplt)
!!find the rest
DO i = 2, N
        CALL NextPoint(follower(i-1),follower(i),intplt)
ENDDO
CALL PGSAVE
CALL PGSCI(2)
CALL PGPT(N,real(follower(:).x),real(follower(:).y),2)
CALL PGUNSA

ENDSUBROUTINE

SUBROUTINE InitStellar
USE FollowMO
USE STELLARDISK
IMPLICIT NONE
type(typspiral)                 ::spiral
DOUBLE PRECISION                ::r,th
DOUBLE PRECISION,POINTER        ::x=>null(),y=>null()
INTEGER                         ::i,j

ALLOCATE(StellarDensity(gas.n,gas.n))
CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)

!$OMP PARALLEL SHARED(StellarDensity) PRIVATE(j,x,y,r,th)
!$OMP DO
DO i = 1, gas.n
DO j = 1, gas.n
        x=>gas.x(i)
        y=>gas.y(j)
        r = sqrt(x**2+y**2)
        th= atan2(y,x)
        StellarDensity(i,j) = sigma1(r,th,spiral)
!       StellarDensity(i,j) = 0.d0
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
CALL spiral.free
ENDSUBROUTINE

SUBROUTINE FreeStellar
USE FollowMO
IMPLICIT NONE
DEALLOCATE(StellarDensity)
ENDSUBROUTINE

SUBROUTINE InitPlot
USE argument
IMPLICIT NONE
INTEGER                                 ::PGOPEN
DOUBLE PRECISION                        ::domain
write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'

IF(xn_set)THEN
        IF(PGOPEN(ACHAR(48+xn)//'/xs').NE.1)STOP
ELSE
        IF (PGOPEN('/xs') .NE. 1) STOP
ENDIF

CALL PGSUBP(-3,2)

ENDSUBROUTINE

SUBROUTINE InitGasPlot
USE FollowMO
USE plotting
USE argument
IMPLICIT NONE
DOUBLE PRECISION                        ::domain
IF(.not.zmax_set)zmax = maxval(gas.density)
domain = maxval(gas.x)
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
CALL meshplot(gas.density,gas.n,domain,zmax,zmin_in=0.d0)

!!plot circle
IF(plot_circ)THEN
        CALL PGSAVE
        CALL PGSCI(0)
        CALL PGSFS(2)
        CALL PGCIRC(0.,0.,3.)
        CALL PGCIRC(0.,0.,4.)
        CALL PGCIRC(0.,0.,5.)
        CALL PGUNSA
ENDIF

ENDSUBROUTINE

SUBROUTINE InitStarPlot
USE FollowMO
USE plotting
USE argument,only:plot_circ
IMPLICIT NONE
INTEGER                                 ::PGOPEN
DOUBLE PRECISION                        ::domain
DOUBLE PRECISION                        ::zmax

zmax = maxval(StellarDensity)
write(6,*)achar(27)//'[33m Plotting z scale of star:',zmax,achar(27)//'[0m'

domain = maxval(gas.x)
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
CALL meshplot(StellarDensity,gas.n,domain,zmax,zmin_in=0.d0)

!!plot circle
IF(plot_circ)THEN
        CALL PGSAVE
        CALL PGSCI(0)
        CALL PGSFS(2)
        CALL PGCIRC(0.,0.,3.)
        CALL PGCIRC(0.,0.,4.)
        CALL PGCIRC(0.,0.,5.)
        CALL PGUNSA
ENDIF

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

DO i = 1, 2
thi = minval(Followers(:,i).th)*0.8
thf = maxval(Followers(:,i).th)*1.1
ampi = minval(Followers(:,i).amp)*0.8
ampf = maxval(Followers(:,i).amp)*1.2
ri  = minval(Followers(:,i).r)*0.8
rf  = maxval(Followers(:,i).r)*1.2

IF(zmax_set)ampf = zmax

CALL PGPANL(2,i)
!!plot th to amp
CALL PGSWIN(real(pi),real(-pi),real(ampi),real(ampf))
CALL PGSVP(0.1,0.9,0.6,0.9)
CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
CALL PGSAVE
CALL PGLINE(N,real(Followers(:,i).th),real(Followers(:,i).amp))
CALL PGUNSA
CALL PGLAB('\gh','Amp','\gh to Amp')

!!plot r to amp
CALL PGSWIN(real(ri),real(rf),real(ampi),real(ampf))
CALL PGSVP(0.1,0.9,0.1,0.4)
CALL PGBOX('BCTSN',0.0,0,'BCTSN',0.0,0)
CALL PGSAVE
CALL PGLINE(N,real(Followers(:,i).r),real(Followers(:,i).amp))
CALL PGUNSA
CALL PGLAB('r','Amp','r to Amp')
ENDDO
ENDSUBROUTINE

SUBROUTINE PlotPhase(N,M,Followers)
USE FollowMO
USE argument,only:output
USE math,only:intplt1
IMPLICIT NONE
CHARACTER(len=225)                    ::CH
TYPE(tyFollower)                     ::Followers(N,M)
INTEGER                              ::N,M
DOUBLE PRECISION                     ::dr,ri,rf,r,dmax,dmin
INTEGER,PARAMETER                    ::P=500
DOUBLE PRECISION                     ::dat(P,2)
INTEGER                              ::i
ri = max(minval(Followers(:,1).r),minval(Followers(:,2).r))
rf = min(maxval(Followers(:,1).r),maxval(Followers(:,2).r))
dr = (rf - ri)/dble(P)

DO i = 1, P
        r = ri + dr*dble(i)
        dat(i,1) = r
        dat(i,2) = intplt1(Followers(:,1).th,Followers(:,1).r,r) &
                 - intplt1(Followers(:,2).th,Followers(:,2).r,r)
ENDDO

dmax = maxval(dat(:,2))*0.8
dmin = minval(dat(:,2))*1.1

CALL PGENV(real(ri)*0.9,real(rf)*1.1, -real(pi)/2.,real(pi)/4.,0,0)
!CALL PGENV(3.5,7.2, -real(pi)/2.,0.,0,0)
CALL PGLINE(P,real(dat(:,1)),real(dat(:,2)))
CALL PGLINE(2,real((/ri,rf/)),real((/-pi/2.d0,-pi/2.d0/)))
CALL PGLAB('Radius (kpc)','Offset (radians)','')


IF(output)THEN
        !!get working path as filename of part of the file name of output.
        CALL GETCWD(CH)
        CH = TRIM(CH(INDEX(CH,'/',.true.)+1:LEN(CH)))
        !!start output
        OPEN(10,file=TRIM(CH)//'_phase.dat')
        write(10,100)'#r','th'
        DO i = 1, P
                write(10,100)dat(i,1),dat(i,2)
        ENDDO
        CLOSE(10)
        100 FORMAT(5(G13.6,8X))
ENDIF

ENDSUBROUTINE

SUBROUTINE FirstPoint(Follower,intplt)
USE FollowMO
IMPLICIT NONE
TYPE(tyFollower)                ::Follower
TYPE(typintplt2)                ::intplt
CHARACTER(1)                    ::CH
DOUBLE PRECISION                ::PGCURS
REAL                            ::rx,ry
INTEGER                         ::ierr

ierr = PGCURS(rx,ry,CH)
CALL Follower.init(dble(rx),dble(ry),intplt)
write(6,'("using : x = ",F5.2,"  y = ",F5.2)')rx,ry

ENDSUBROUTINE

SUBROUTINE NextPoint(Follower,NxFollower,intplt)
USE FollowMO
USE math,only:typintplt2
IMPLICIT NONE
TYPE(tyFollower),INTENT(IN)     ::Follower
TYPE(tyFollower),INTENT(OUT)    ::NxFollower
TYPE(typintplt2)                ::intplt
DOUBLE PRECISION,PARAMETER      ::r=0.1d0
DOUBLE PRECISION                ::th,dth,amp
INTEGER,PARAMETER               ::N = 200
DOUBLE PRECISION                ::dat(N,2)
INTEGER                         ::i


dth = -pi/dble(N) 
DO i = 1, N
        th = dth*dble(i)+ atan2(Follower.y,Follower.x)
        amp= intplt.find(Follower.x+r*cos(th),Follower.y+r*sin(th))
        dat(i,1) = th
        dat(i,2) = amp
ENDDO
th = dat(maxloc(dat(:,2),1),1)
CALL NxFollower.init(Follower.x+r*cos(th),Follower.y+r*sin(th),intplt)
ENDSUBROUTINE

SUBROUTINE OutputPhase(N,M,Followers)
USE FollowMO
IMPLICIT NONE
TYPE(tyFollower)                       ::Followers(N,M)
INTEGER                                ::N,M
INTEGER                                ::PGOPEN

IF(PGOPEN('PhaseShift.ps/ps').EQ.0)STOP
CALL PlotPhase(N,M,Followers)
ENDSUBROUTINE
