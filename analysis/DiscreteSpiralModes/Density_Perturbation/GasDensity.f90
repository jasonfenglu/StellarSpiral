module projections
        DOUBLE PRECISION,PARAMETER      ::pi_n = atan(1.d0)*4.d0
        DOUBLE PRECISION,SAVE           ::aline,argaline
        DOUBLE PRECISION,SAVE           ::apitc
        DOUBLE PRECISION,SAVE           ::aolin
CONTAINS 
FUNCTION R(th)
IMPLICIT NONE
DOUBLE PRECISION                ::R(2,2)
DOUBLE PRECISION                ::th
R(1,1) = cos(th)
R(1,2) = -sin(th)
R(2,1) = sin(th)
R(2,2) = cos(th)
ENDFUNCTION

FUNCTION C(th)
IMPLICIT NONE
DOUBLE PRECISION                ::C(2,2)
DOUBLE PRECISION                ::th
C      = 0.d0
C(1,1) = cos(th)
C(2,2) = 1.d0
ENDFUNCTION

FUNCTION INVC(th)
IMPLICIT NONE
DOUBLE PRECISION                ::INVC(2,2)
DOUBLE PRECISION                ::th
INVC      = 0.d0
INVC(1,1) = 1.d0/cos(th)
INVC(2,2) = 1.d0
ENDFUNCTION
endmodule

MODULE argument
IMPLICIT  NONE
        DOUBLE PRECISION        ::zmax      = 30.d0
        CHARACTER(len=225)      ::filename
        INTEGER                 ::frames(2)
        INTEGER                 ::drawmode  = 1         !1 for draw static
                                                        !density, 2 for movie,
        LOGICAL                 ::toproject = .false.
        LOGICAL                 ::rauto     = .true.
        LOGICAL                 ::drawcir   = .false.
        LOGICAL                 ::drawstellar=.false.
        LOGICAL                 ::givefnm   = .false.
        LOGICAL                 ::drawq     = .false.
        LOGICAL                 ::save      = .false.
        LOGICAL                 ::cut       = .false.
CONTAINS
SUBROUTINE readarg
USE projections,only:argaline
IMPLICIT NONE
CHARACTER(len=32)               ::arg
INTEGER                         ::i
INTEGER                         ::ioerr
!read in rotating parameters and gas filename by stdin
if(iargc().ne.0)then
        DO i = 1, iargc()
                CALL getarg(i,arg)
                SELECT CASE(arg)
                CASE('--help','-h')
                        write(6,'(a)')'Draw gas density distribution.       '
                        write(6,'(a)')'usage:   GasDensity.exe [option]            '
                        write(6,'(a)')'options:                         '
                        write(6,'(a)')'         -c, --circle            To draw circles.       '
                        write(6,'(a)')'         -i [file name.h5]       Read file.'
                        write(6,'(a)')'         -m [start] [end]        Read files.'
                        write(6,'(a)')'         -p, --project [degree]  To project the density.'
                        write(6,'(a)')'         -r, --zauto             Automatic find scale of z'
                        write(6,'(a)')'         -z, --zmax [scale of z] Specified the scale of z'
                        write(6,'(a)')'         -h, --help              Show this help page     '
                        write(6,'(a)')'         -s, --stellar           Draw stellar contour.   '
                        write(6,'(a)')'         -q,                     Draw instability.'
                        write(6,'(a)')'         -v,                     Version information.'
                        write(6,'(a)')'         --save                  Save as png file.   '
                        write(6,'(a)')'         --cut                   Make a cut at CO.   '
                        STOP
                CASE('--circle','-c')
                        drawcir = .true.
                CASE('-i')
                        CALL getarg(i+1,arg)
                        READ(arg,*)filename
                        givefnm = .true.
                CASE('-m')
                        CALL getarg(i+1,arg)
                        READ(arg,*)frames(1)
                        CALL getarg(i+2,arg)
                        READ(arg,*)frames(2)
                        drawmode = 2
                CASE('--project','-p')
                        toproject = .true.
                        CALL getarg(i+1,arg)
                        READ(arg,*,iostat=ioerr)argaline
                        if(ioerr.ne.0)argaline = -3.d0
                CASE('--stellar','-s')
                        drawstellar = .true.
                CASE('--save')
                        save = .true.
                CASE('--zmax','-z')
                        CALL getarg(i+1,arg)
                        READ(arg,*)zmax
                        rauto = .false.
                CASE('-q')
                        drawq = .true.
                CASE('--cut')
                        cut = .true.
                CASE('-v')
                        write(6,'(a)')'Compiled at: '//__DATE__//' '//__TIME__
                        STOP
                ENDSELECT
        ENDDO
ENDIF
ENDSUBROUTINE
ENDMODULE

PROGRAM density1
USE PLOTTING
USE STELLARDISK_MODEL
USE STELLARDISK,only:FindSpiral,pi_n=>pi,sigma1
USE projections,only:argaline
USE io
USE argument
IMPLICIT NONE
TYPE gastyp
        CHARACTER(len=32)           ::filename
        DOUBLE PRECISION,ALLOCATABLE::density(:,:)
        DOUBLE PRECISION,ALLOCATABLE::x(:)
        DOUBLE PRECISION,ALLOCATABLE::y(:)
        DOUBLE PRECISION            ::dx,dy
ENDTYPE
type(gastyp)                    ::gas
INTEGER                         ::i,j,k,l
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::domain= 10.d0,dx,dy,r,th,pf(2),pi(2)
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::GasDensity(:,:)
DOUBLE PRECISION                ::limit = 100.d0
DOUBLE PRECISION                ::d
DOUBLE PRECISION                ::intb(4)
INTEGER,PARAMETER               ::n=600
INTEGER                         ::PGBEG
type(typspiral)                 ::spiral


!!Reading from arguments
CALL readarg
gas.filename = filename

!!check if gas reading files are set:
if(.not.(givefnm.or.drawmode.eq.2))then
        write(0,'(a)')'No input file or files set'
        STOP
ENDIF

!preparing spiral
CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)

!set up grid
ALLOCATE(density(n,n,3))!1 for stellar, 2 for gas, 3 for q of gas
ALLOCATE(xcoord(n))
ALLOCATE(ycoord(n))

dx = domain/dble(n)*2.d0
dy = domain/dble(n)*2.d0
DO i = 1, n
        xcoord(i) =  dx*0.5d0 - domain + dble(i-1)*dx
        ycoord(i) =  dy*0.5d0 - domain + dble(i-1)*dy
ENDDO

!!Filling stellar density
!$OMP PARALLEL SHARED(density,spiral,gas) PRIVATE(j,r,th,pi,pf,d,k,l,intb)
!$OMP DO 
DO i = 1, n
DO j = 1, n
        !read coordinate
        pf = (/xcoord(i),ycoord(j)/)
        pi = pf
        !project to original coordinate if nessesary
        if(toproject)CALL projection(pi,pf)
        !Fill in Stellar Density
        r  = sqrt(pi(1)**2+pi(2)**2)
        th = atan2(pi(2),pi(1))
        IF(r.lt.spiral.fortoone)then
                d  = sigma1(r,th,spiral)
        ELSE    
                d  = 0.d0
        ENDIF
        !remove wired pixel
        if(isnan(d))d = 0.d0
        density(i,j,1) = d
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 
FORALL(i = 1:n,j = 1:n,density(i,j,1)<0)
        density(i,j,1) = 0.d0
ENDFORALL


points(1,:) = (/0.0,10.636/)
points(2,:) = (/0.0,-10.636/)
points(3,:) = (/0.0,4.727/)
points(4,:) = (/0.0,-4.727/)
CALL dprojection(points)

!!looping to draw movie
!read in gas density
SELECT CASE(drawmode)
CASE(1) ! plot static density map
        CALL ReadGasDensity
        CALL FillGasDensity
        IF(.not.drawq)THEN              !plot density only
                IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
                !!set page
                CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
                CALL PlotGasDensity(gas.filename)
        ELSE                            !plot density and q
                IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
                !!set coordinate
                CALL PGSUBP(2,1)
                CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
                CALL PlotGasDensity(gas.filename)
                CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
                CALL Plotq(gas.filename)
        ENDIF
CASE(2) ! plot movie
        IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
        !!set plot layout
        IF(.not.drawq)THEN
                CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        ELSE
                CALL PGSUBP(2,1)
                CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        ENDIF
        DO i = frames(1),frames(2)
        write(gas.filename,'("M",I0.4,".h5")')i
        print *,gas.filename
        CALL ReadGasDensity
        CALL FillGasDensity
        IF(rauto)THEN
                CALL meshplot(density(:,:,2),n,domain,maxval(density(:,:,2)),zmin_in=0.d0)
        ELSE
                CALL meshplot(density(:,:,2),n,domain,zmax,zmin_in=0.d0)
        ENDIF
        CALL PGETXT
        CALL PGPTXT(5.,8.,0.,0.,gas.filename)
        CALL PlotCircle
        CALL PGBBUF
        CALL PlotStellarDensity(0)
        CALL PGEBUF
        IF(drawq)THEN
                CALL PGPANL(2,1)
                CALL Plotq(gas.filename)
                CALL PGPANL(1,1)
        ELSE

        ENDIF
        ENDDO
ENDSELECT

CALL PGCLOS

CALL PrintGasDensity
!CALL StellarGasPlot(density,n,domain,4)
!!close h5 files writing
!CALL h5write(xcoord,N,'density.h5','xcoord')
!CALL h5write(ycoord,N,'density.h5','ycoord')
!CALL h5write(density(:,:,2),N,N,'density.h5','density')

DEALLOCATE(xcoord)
DEALLOCATE(ycoord)
DEALLOCATE(density)
!CALL PhaseIntegrate
!CALL ENDSTELLARDISK
DEALLOCATE(gas.density)
DEALLOCATE(gas.x)
DEALLOCATE(gas.y)
STOP

CONTAINS

SUBROUTINE ReadGasDensity
IMPLICIT NONE
INTEGER                         ::fsize(2)

fsize = h5size2d(gas.filename,'density')
IF(.not.ALLOCATED(gas.density))ALLOCATE(gas.density(fsize(1),fsize(2)))
IF(.not.ALLOCATED(gas.x))ALLOCATE(gas.x(fsize(1)))
IF(.not.ALLOCATED(gas.y))ALLOCATE(gas.y(fsize(2)))

CALL h5read(gas.density,fsize(1),fsize(2),gas.filename,'density')
CALL h5read(gas.x,fsize(1),gas.filename,'x')
CALL h5read(gas.y,fsize(2),gas.filename,'y')
gas.dx = gas.x(2)-gas.x(1)
gas.dy = gas.y(2)-gas.y(1)

ENDSUBROUTINE

SUBROUTINE PlotGasDensity(title)
USE argument
USE plotting
IMPLICIT NONE
CHARACTER(len=32)                       ::title
INTEGER                                 ::i

IF(rauto)THEN
        zmax = maxval(density(:,:,2))*1.1d0
        write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'
ELSE
        write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'
ENDIF

!!plot gas 
CALL meshplot(density(:,:,2),n,domain,zmax,zmin_in=0.d0)
!!plot circle
CALL PlotCircle

!!plot stellar contour
CALL PlotStellarDensity(0)
        
CALL PGLAB('kpc','kpc',title)

ENDSUBROUTINE

SUBROUTINE PrintGasDensity
USE argument
USE plotting
IMPLICIT NONE
INTEGER                                 ::i
INTEGER                                 ::PGOPEN
CHARACTER(len=225)                      ::CH


CALL GETCWD(CH)
CH = TRIM(CH(INDEX(CH,'/',.true.)+1:LEN(CH)))

IF (PGOPEN(TRIM(CH)//'_GasDensity.png/png') .NE. 1) STOP
!!plot gas 
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
CALL meshplot(density(:,:,2),n,domain,zmax,zmin_in=0.d0)
!!plot circle
CALL PlotCircle

!!plot stellar contour
CALL PlotStellarDensity(0)

CALL PGLAB('kpc','kpc','Gas Density')
CALL PGCLOS

IF(drawq)THEN
        IF (PGOPEN(TRIM(CH)//'_GasDensity_Q.png/png') .NE. 1) STOP
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        CALL plotq(CH)
        CALL PGCLOS
ENDIF

ENDSUBROUTINE

SUBROUTINE plotq(title)
USE argument
USE plotting
IMPLICIT NONE
CHARACTER(*)                            ::title
INTEGER                                 ::i

!IF(rauto)THEN
!        zmax = maxval(density(:,:,3))*1.1d0
!        write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'
!ELSE
!        write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'
!ENDIF

!!plot q of gas 
CALL meshplot(density(:,:,3),n,domain, 3.d0,zmin_in=0.d0)

!!plot circle
CALL PlotCircle

!!plot stellar contour
CALL PlotStellarDensity(1)
        
CALL PGLAB('kpc','kpc','-log(Q)')

ENDSUBROUTINE

SUBROUTINE FillGasDensity
USE STELLARDISK, only:kappa,GravConst,pi_n=>pi
USE galaxy, only:gasdensity
IMPLICIT NONE
!!Filling stellar density
!$OMP PARALLEL SHARED(density,spiral,gas) PRIVATE(j,r,th,pi,pf,d,k,l,intb)
!$OMP DO 
DO i = 1, n
DO j = 1, n
        !read coordinate
        pf = (/xcoord(i),ycoord(j)/)
        pi = pf
        !project to original coordinate if nessesary
        if(toproject)CALL projection(pi,pf)
        r  = sqrt(pi(1)**2+pi(2)**2)
        th = atan2(pi(2),pi(1))
        !Fill in Gas Density
        !find interpolating grid
        !!if pixel is not in simulation grid
        IF(pi(1)<minval(gas.x).or.pi(1)>maxval(gas.x).or.pi(2)<minval(gas.y).or.pi(2)>maxval(gas.y))then
                density(i,j,2) = 0.d0
        ELSE
                r = sqrt(pi(1)**2+pi(2)**2)
                DO k = 1, size(gas.x)
                        if(gas.x(k)>pi(1))EXIT
                ENDDO
                DO l = 1, size(gas.y)
                        if(gas.y(l)>pi(2))EXIT
                ENDDO
                !change k,l to the index for smaller neibor 
                k = k - 1
                l = l - 1
                !remap pi to unit square
                pi(1) = (pi(1)-gas.x(k))/gas.dx
                pi(2) = (pi(2)-gas.y(l))/gas.dy
                intb(1) = gas.density(k,l)
                intb(2) = gas.density(k+1,l) - intb(1)
                intb(3) = gas.density(k,l+1) - intb(1)
                intb(4) = 0.d0
                intb(4) = -sum(intb(:)) + gas.density(k+1,l+1)
                density(i,j,2) = intb(1) + intb(2)*pi(1) + intb(3)*pi(2) &
                               + intb(4)*pi(1)*pi(2)
                density(i,j,3) = &
                kappa(r,spiral)*8.d0/pi_n/GravConst/(density(i,j,2))/1d6
                density(i,j,3) = -log(density(i,j,3))
                
        ENDIF
        IF(cut)THEN
                IF(r>spiral.co)THEN
                        density(i,j,2) = 0.d0
                        density(i,j,3) = 0.d0
                ENDIF
        ENDIF
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 

ENDSUBROUTINE

SUBROUTINE PlotStellarDensity(color)
USE argument,only:drawstellar
USE plotting
IMPLICIT NONE
INTEGER                                 ::color
IF(drawstellar)then
        write(6,*)achar(27)//'[33m Drawing stellar contour',achar(27)//'[0m'
        !CALL PGSCI(0)
        CALL PGSAVE
        CALL PGSCI(color)
        CALL contourplot(density(:,:,1),n,domain,4)
        CALL PGUNSA
ENDIF

ENDSUBROUTINE

END PROGRAM

SUBROUTINE projection(pi,pf)
USE projections
IMPLICIT NONE
DOUBLE PRECISION                ::pi(2),pf(2)
!!BLAS
CHARACTER(1)                    ::TRANS
DOUBLE PRECISION                ::ALPHA = 1.d0
DOUBLE PRECISION                ::A(2,2)
DOUBLE PRECISION                ::X(2),Y(2)
DOUBLE PRECISION                ::BETA  = 0.d0
INTEGER                         ::M = 2
INTEGER                         ::N = 2
INTEGER                         ::LDA = 2
INTEGER                         ::INCX = 1
INTEGER                         ::INCY = 1

CALL set_angles

X = pf
A = R(-aolin)
TRANS  = 'n'
CALL DGEMV (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
X = Y
A = INVC(apitc)
TRANS  = 'n'
CALL DGEMV (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
X = Y
A = R(-aline)
TRANS  = 'n'
CALL DGEMV (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)

pi = y
ENDSUBROUTINE

SUBROUTINE dprojection(p)
USE projections
USE plotting,only:points
IMPLICIT NONE
DOUBLE PRECISION                ::pi(2),pf(2)
REAL                            ::p(4,2)
INTEGER                         ::i
!!BLAS
CHARACTER(1)                    ::TRANS
DOUBLE PRECISION                ::ALPHA = 1.d0
DOUBLE PRECISION                ::A(2,2)
DOUBLE PRECISION                ::X(2),Y(2)
DOUBLE PRECISION                ::BETA  = 0.d0
INTEGER                         ::M = 2
INTEGER                         ::N = 2
INTEGER                         ::LDA = 2
INTEGER                         ::INCX = 1
INTEGER                         ::INCY = 1

call set_angles

DO i = 1,4
        X = points(i,:)
        A = R(aolin)
        TRANS  = 'n'
        CALL DGEMV (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
        p(i,:) = y
ENDDO



ENDSUBROUTINE

SUBROUTINE set_angles
USE projections
IMPLICIT NONE
!observed angle of line of node 28 degree
 aolin = -28.3d0/180.d0*pi_n
!aolin = -32.d0/180.d0*pi_n
!inclination angle is 55 degree ?
apitc = 55.d0/180.d0*pi_n
!tuned angle of line of node
!aline = 30.d0/180.d0*pi_n
aline =  argaline/180.d0*pi_n
ENDSUBROUTINE

SUBROUTINE PhaseIntegrate(spiral)
USE STELLARDISK_MODEL
USE STELLARDISK,only:k3sqrt,pi
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION                ::ri,rf
DOUBLE PRECISION                ::RE,AE,RR,ERR,ANS
INTEGER                         ::IFLAG
ri = 0.d0
rf = 4.d0
RR = 2.d0
RE = 1d-8
AE = 1d-8
CALL DFZERO(F,ri,rf,RR,RE,AE,IFLAG)
ri = 0.d0
rf = min(5.1d0,rf)
ERR = 1.d-8
CALL DGAUS8(g,ri,rf,ERR,ANS,IFLAG)
print *,'phase before zero',ans/pi,'to',rf
CONTAINS
FUNCTION F(r)
IMPLICIT NONE
DOUBLE PRECISION                ::r,F
F = k3sqrt(r,spiral)
ENDFUNCTION
FUNCTION g(r)
IMPLICIT NONE
DOUBLE PRECISION                ::r,G
g = REAL(sqrt(k3sqrt(r,spiral)))
ENDFUNCTION
ENDSUBROUTINE

SUBROUTINE SOLVE(rho,phi,N,L)
IMPLICIT NONE
DOUBLE PRECISION        ::rho(2*n,2*n),phi(2*n,2*n),L
REAL,ALLOCATABLE        ::rrho(:,:)
real        ::A,B,C,D,ELMBDA,PERTRB
real,ALLOCATABLE::BDA(:),BDB(:),BDC(:),BDD(:),w(:)
INTEGER                 ::IDIMF,MBDCND,NBDCND
INTEGER                 ::IERROR,wdim,n,i


IDIMF = 2*n
A = -real(l)
B = real(l)
MBDCND = 1
C = -real(l)
D = real(l)
NBDCND = 1
ELMBDA = 0.0
wdim = 13 + INT(LOG(dble(N))/log(2.))*N + 4*N
wdim = wdim*4

ALLOCATE(rrho(2*n,2*n))
ALLOCATE(BDA(2*n))
ALLOCATE(BDB(2*n))
ALLOCATE(BDC(2*n))
ALLOCATE(BDD(2*n))
ALLOCATE(W(wdim))

!print *,n,size(rho),size(rrho)
rrho = real(rho)
BDA(:) = 0.
BDB = BDA
BDC = BDA
BDD = BDA

CALL HSTCRT (A, B, 2*N, MBDCND, BDA, BDB, C, D, 2*N, NBDCND, BDC, BDD, ELMBDA, rrho,IDIMF, PERTRB, IERROR, W)
phi = dble(rrho)

DEALLOCATE(rrho)
DEALLOCATE(BDA)
DEALLOCATE(BDB)
DEALLOCATE(BDC)
DEALLOCATE(BDD)
DEALLOCATE(W)

ENDSUBROUTINE

SUBROUTINE PlotCircle
USE argument
IMPLICIT NONE
IF(.not.toproject.and.drawcir)THEN
        write(6,*)achar(27)//'[33m Drawing circles',achar(27)//'[0m'
        CALL PGSAVE
        CALL PGSFS(2)
        CALL PGSCI(0)
        CALL PGCIRC(0.,0.,1.26)
        CALL PGCIRC(0.,0.,2.36)
        CALL PGCIRC(0.,0.,4.72)
        CALL PGCIRC(0.,0.,10.636)
        CALL PGCIRC(0.,0.,8.83)
        CALL PGUNSA
ENDIF
ENDSUBROUTINE

