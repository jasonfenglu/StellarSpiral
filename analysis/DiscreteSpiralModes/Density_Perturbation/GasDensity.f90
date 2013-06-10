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
        INTEGER                 ::frames(2)
        LOGICAL                 ::toproject = .false.
        LOGICAL                 ::rauto     = .true.
        LOGICAL                 ::drawcir   = .false.
        LOGICAL                 ::drawstellar=.false.
        LOGICAL                 ::movie     = .false.
        LOGICAL                 ::givefnm   = .false.
ENDMODULE

PROGRAM density1
USE PLOTTING
USE STELLARDISK_MODEL
USE STELLARDISK,only:FindSpiral,pi_n=>pi,sigma1,phi1,FindPhi1,SpiralForce
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
INTEGER                         ::ioerr
type(typspiral)                 ::spiral

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
                        STOP
                CASE('--circle','-c')
                        drawcir = .true.
                CASE('-i')
                        CALL getarg(i+1,arg)
                        READ(arg,*)gas.filename
                        givefnm = .true.
                CASE('-m')
                        CALL getarg(i+1,arg)
                        READ(arg,*)frames(1)
                        CALL getarg(i+2,arg)
                        READ(arg,*)frames(2)
                        movie = .true.
                CASE('--project','-p')
                        toproject = .true.
                        CALL getarg(i+1,arg)
                        READ(arg,*,iostat=ioerr)argaline
                        if(ioerr.ne.0)argaline = -3.d0
                CASE('--stellar','-s')
                        drawstellar = .true.
                CASE('--zmax','-z')
                        CALL getarg(i+1,arg)
                        READ(arg,*)zmax
                        rauto = .false.
                ENDSELECT
        ENDDO
ENDIF

!!check if gas reading files are set:
if(.not.(givefnm.or.movie))then
        write(0,'(a)')'No input file or files set'
        STOP
ENDIF

!preparing spiral
CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)


!set up grid
ALLOCATE(density(n,n,2))
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


points(1,:) = (/0.0,10.636/)
points(2,:) = (/0.0,-10.636/)
points(3,:) = (/0.0,4.727/)
points(4,:) = (/0.0,-4.727/)
CALL dprojection(points)

!!looping to draw movie
!read in gas density
IF(movie)THEN
        DO i = frames(1),frames(2)
                write(gas.filename,'("M",I0.4,".h5")')i
                print *,gas.filename
        CALL ReadGasDensity
        CALL FillGasDensity
        FORALL(i = 1:n,j = 1:n,density(i,j,1)<0)
                density(i,j,1) = 0.d0
        ENDFORALL
        CALL plot(density,n,domain,gas.filename)
        ENDDO
ELSE
        CALL ReadGasDensity
        CALL FillGasDensity
        FORALL(i = 1:n,j = 1:n,density(i,j,1)<0)
                density(i,j,1) = 0.d0
        ENDFORALL

        !!Check if pspd is true
        print *,'error',abs(spiral.error)
        if(abs(spiral.error).gt.1d-5)then
                write(0,*)'!!!!!! wrong pspd:'
                write(0,*)'error:',abs(spiral.error)
                write(0,*)'pspd:',spiral.w
        endif
        CALL plot(density,n,domain,gas.filename)
ENDIF

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

SUBROUTINE plot(F,n,domain,title)
USE argument
USE plotting
IMPLICIT NONE
CHARACTER(len=32)                       ::title
DOUBLE PRECISION,INTENT(IN)             ::F(:,:,:)      !plotting data
DOUBLE PRECISION,INTENT(IN)             ::domain        !plot range
INTEGER,INTENT(IN)                      ::n             !dimentsion
INTEGER                                 ::PGBEG
INTEGER                                 ::i

IF(rauto)THEN
        zmax = maxval(F(:,:,2))*1.1d0
        write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'
ELSE
        write(6,*)achar(27)//'[33m Plotting z scale :',zmax,achar(27)//'[0m'
ENDIF

IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP

!!plot on screen and save as png file
DO i = 1, 2
!!set coordinate
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
!!plot gas 
CALL meshplot(F(:,:,2),n,domain,zmax,zmin_in=0.d0)
!!plot circle
IF(.not.toproject.and.drawcir)THEN
        write(6,*)achar(27)//'[33m Drawing circles',achar(27)//'[0m'
        CALL PGSFS(2)
        CALL PGSCI(0)
        CALL PGCIRC(0.,0.,1.26)
        CALL PGCIRC(0.,0.,2.36)
        CALL PGCIRC(0.,0.,4.72)
        CALL PGCIRC(0.,0.,10.636)
        CALL PGCIRC(0.,0.,8.83)
ENDIF

!!plot stellar contour
IF(drawstellar)then
        write(6,*)achar(27)//'[33m Drawing stellar contour',achar(27)//'[0m'
        CALL PGSCI(0)
        CALL contourplot(F(:,:,1),n,domain,4)
        CALL PGSCI(1)
ENDIF
        
CALL PGLAB('kpc','kpc',title)

CALL PGCLOS
IF (PGBEG(0,'density.png/png',1,1) .NE. 1) STOP
ENDDO

ENDSUBROUTINE

SUBROUTINE FillGasDensity
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
                
        ENDIF
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 

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

