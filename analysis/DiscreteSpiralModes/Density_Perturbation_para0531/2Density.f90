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

PROGRAM caldensity2
USE PLOTTING
USE STELLARDISK_MODEL
USE STELLARDISK,only:FindSpiral,pi_n=>pi,sigma1,k3sqrtlog
USE projections,only:argaline
IMPLICIT NONE
INTEGER                         ::i,j,k
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::domain= 12.d0,dx,dy,r,th,pf(2),pi(2)
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::potential(:,:)
DOUBLE PRECISION,ALLOCATABLE    ::force(:,:,:)
DOUBLE PRECISION,ALLOCATABLE    ::k3(:,:),u(:,:),h(:,:)
DOUBLE PRECISION                ::limit = 100.d0
DOUBLE PRECISION                ::d
DOUBLE PRECISION                ::modeth(2)
INTEGER,PARAMETER               ::n=800
INTEGER                         ::nmode
type(typspiral)                 ::spiral(2)
LOGICAL                         ::toproject
namelist /densitypara/ toproject

open(10,file='para.list')
read(10,nml=densitypara)
close(10)

SELECT CASE(iargc())
CASE(1)
        CALL getarg(1,arg)
        READ(arg,*)argaline
        print *,'read in argaline for n0 = ',argaline
CASE(2)
        CALL getarg(1,arg)
        READ(arg,*)argaline
        print *,'read in argaline for n0 = ',argaline
        CALL getarg(2,arg)
        READ(arg,*)modeth(2)
        print *,'read in argaline for n1 = ',modeth(2)
        
CASE DEFAULT 
        argaline = 3.d0
ENDSELECT

!Setup grid
dx = domain/dble(n)
dy = domain/dble(n)

ALLOCATE(density(2*n,2*n,2))
ALLOCATE(xcoord(2*n))
ALLOCATE(ycoord(2*n))

DO i = 1, n*2
        xcoord(i) =  dx*0.5d0 - domain + dble(i)*dx
        ycoord(i) =  dy*0.5d0 - domain + dble(i)*dx
ENDDO

CALL stdpara.readstd
DO i = 1, 2
        CALL spiral(i).init(500,12.d0,stdpara,i)
        CALL spiral(i).readw(i)
        CALL FindSpiral(spiral(i))
        CALL k3sqrtlog(spiral(i))
ENDDO

modeth(1) = 0.d0
modeth(2) = 55.d0/180.d0*pi_n

!mode loop
DO nmode = 1, 2
!$OMP PARALLEL SHARED(density) PRIVATE(j,r,th,pi,pf,d) FIRSTPRIVATE(spiral)
!$OMP DO 
DO i = 1, n*2
DO j = 1, n*2
        pf = (/xcoord(i),ycoord(j)/)
        pi = pf
        if(toproject)CALL projection(pi,pf)
        r  = sqrt(pi(1)**2+pi(2)**2)
        th = atan2(pi(2),pi(1)) + modeth(nmode)
        d  = sigma1(r,th,spiral(nmode))
        !===================
        !ignore d too high
        if(r.gt.12.d0)d = 0.d0
        !ignore inside r=1.26  
!       if(r.lt.2.26)d = d*(1.d0 - cos(r/2.26d0*pi_n/2.d0))
        !ignore d that is not exist during coordinate transformation
        if(isnan(d))d = 0.d0
        !ignore value below detection limit
!       if(abs(d).lt.limit)d = 0.d0
        density(i,j,nmode) = d
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 
ENDDO

if(maxval(abs(spiral(:).error)).gt.1d-5)write(0,*)'!!!!!! wrong pspd'

points(1,:) = (/0.0,10.636/)
points(2,:) = (/0.0,-10.636/)
points(3,:) = (/0.0,4.727/)
points(4,:) = (/0.0,-4.727/)
CALL dprojection(points)
ALLOCATE(k3(spiral(1).n,2))
k3(:,1) = real(spiral(1).k3)
k3(:,2) = real(spiral(2).k3)
ALLOCATE(u(spiral(1).n,2))
u(:,1) = abs(spiral(1).u(2,:))
u(:,2) = abs(spiral(2).u(2,:))
ALLOCATE(h(spiral(1).n,2))
h(:,1) = abs(spiral(1).h1(:))
h(:,2) = abs(spiral(2).h1(:))
CALL plotdensity2(density,n,domain,spiral(1).r,k3,spiral(1).n,u)
DEALLOCATE(k3)
DEALLOCATE(u)
DEALLOCATE(h)

!DEALLOCATE(potential)
DEALLOCATE(xcoord)
DEALLOCATE(ycoord)
DEALLOCATE(density)
CALL PhaseIntegrate
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

SUBROUTINE PhaseIntegrate
!USE STELLARDISK,only:k3sqrt,pi
!IMPLICIT NONE
!DOUBLE PRECISION                ::ri,rf
!DOUBLE PRECISION                ::RE,AE,RR,ERR,ANS
!INTEGER                         ::IFLAG
!ri = 0.d0
!rf = 4.d0
!RR = 2.d0
!RE = 1d-8
!AE = 1d-8
!CALL DFZERO(F,ri,rf,RR,RE,AE,IFLAG)
!ri = 0.d0
!rf = min(5.1d0,rf)
!ERR = 1.d-8
!CALL DGAUS8(g,ri,rf,ERR,ANS,IFLAG)
!print *,'phase before zero',ans/pi,'to',rf
!CONTAINS
!FUNCTION F(r)
!IMPLICIT NONE
!DOUBLE PRECISION                ::r,F
!F = k3sqrt(r)
!ENDFUNCTION
!FUNCTION g(r)
!IMPLICIT NONE
!DOUBLE PRECISION                ::r,G
!g = REAL(sqrt(k3sqrt(r)))
!ENDFUNCTION
ENDSUBROUTINE
