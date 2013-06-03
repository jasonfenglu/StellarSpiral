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

PROGRAM density1
USE PLOTTING
USE STELLARDISK_MODEL
USE STELLARDISK,only:FindSpiral,pi_n=>pi,sigma1,phi1,FindPhi1,SpiralForce
USE projections,only:argaline
USE io
IMPLICIT NONE
INTEGER                         ::i,j,k
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::domain= 10.d0,dx,dy,r,th,pf(2),pi(2)
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::potential(:,:),potential2(:,:)
DOUBLE PRECISION,ALLOCATABLE    ::force(:,:,:)
DOUBLE PRECISION                ::limit = 100.d0
DOUBLE PRECISION                ::d
INTEGER,PARAMETER               ::n=512
type(typspiral)                 ::spiral
LOGICAL                         ::toproject
namelist /densitypara/ toproject

!read in density plot related options
open(10,file='para.list')
read(10,nml=densitypara)
close(10)

!rotate with input parameters
if(iargc().eq.1)then
        CALL getarg(1,arg)
        READ(arg,*)argaline
        print *,'read in argaline = ',argaline
        else 
        argaline = 3.d0
endif


CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)
dx = domain/dble(n)
dy = domain/dble(n)

ALLOCATE(density(2*n,2*n))
ALLOCATE(xcoord(2*n))
ALLOCATE(ycoord(2*n))

DO i = 1, n*2
        xcoord(i) =  dx*0.5d0 - domain + dble(i)*dx
        ycoord(i) =  dy*0.5d0 - domain + dble(i)*dx
ENDDO

!$OMP PARALLEL SHARED(density,spiral) PRIVATE(j,r,th,pi,pf,d)
!$OMP DO 
DO i = 1, n*2
DO j = 1, n*2
        pf = (/xcoord(i),ycoord(j)/)
        pi = pf
        if(toproject)CALL projection(pi,pf)
        r  = sqrt(pi(1)**2+pi(2)**2)
        th = atan2(pi(2),pi(1))
        IF(r.lt.spiral.fortoone)then
                d  = sigma1(r,th,spiral)
        ELSE    
                d  = 0.d0
        ENDIF
        !===================
        !ignore d too high
        !if(r.gt.12.d0)d = 0.d0
        !ignore d that is not exist during coordinate transformation
        if(isnan(d))d = 0.d0
        density(i,j) = d
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 

!!Check if pspd is true
print *,'error',abs(spiral.error)
if(abs(spiral.error).gt.1d-5)then
        write(0,*)'!!!!!! wrong pspd:'
        write(0,*)'error:',abs(spiral.error)
        write(0,*)'pspd:',spiral.w
endif

if(.false.)then
!!!Find 2d Potential, not using.
ALLOCATE(potential(2*n,2*n))
ALLOCATE(force(n*2,n*2,2))
!$OMP PARALLEL SHARED(force,potential,spiral) PRIVATE(j,r,th,pi,pf,d)
!$OMP DO 
DO i = 1, n*2
DO j = 1, n*2
        pf = (/xcoord(i),ycoord(j)/)
        pi = pf
        if(toproject)CALL projection(pi,pf)
        r  = sqrt(pi(1)**2+pi(2)**2)
        th = atan2(pi(2),pi(1))
        d  = phi1(r,th,spiral,0.d0)
!       force(i,j,:) = SpiralForce(r,th,spiral)
        !===================
        !ignore d too high
        if(r.gt.12.d0)d = 0.d0
        !ignore d that is not exist during coordinate transformation
        if(isnan(d))d = 0.d0
!       potential(i,j) = d
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 
endif
call spiral.printu
call spiral.printh1
points(1,:) = (/0.0,10.636/)
points(2,:) = (/0.0,-10.636/)
points(3,:) = (/0.0,4.727/)
points(4,:) = (/0.0,-4.727/)
CALL dprojection(points)

FORALL(i = 1:n*2,j = 1:n*2,density(i,j)<0)
        density(i,j) = 0.d0
ENDFORALL
CALL countour(density,n*2,domain,4)
!CALL h5io(xcoord,2*N,'density.h5','xcoord')
!CALL h5io(ycoord,2*N,'density.h5','ycoord')
!CALL h5io(density,2*N,2*N,'density.h5','density')

DEALLOCATE(xcoord)
DEALLOCATE(ycoord)
DEALLOCATE(density)
!CALL PhaseIntegrate
!CALL ENDSTELLARDISK
STOP


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
!pitch angle is 55 degree ?
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
