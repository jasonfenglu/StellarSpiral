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

PROGRAM caldensity
USE PLOTTING
USE STELLARDISK,pi_n=>pi
USE projections,only:argaline
IMPLICIT NONE
INTEGER                         ::i,j,k
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::domain= 12.d0,dx,dy,r,th,pf(2),pi(2)
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::potential(:,:)
DOUBLE PRECISION,ALLOCATABLE    ::force(:,:,:)
DOUBLE PRECISION                ::limit = 100.d0
DOUBLE PRECISION                ::d
INTEGER,PARAMETER               ::n=500
type(spiral_type)               ::shared_spiral
LOGICAL                         ::toproject
namelist /densitypara/ toproject

open(10,file='para.list')
read(10,nml=densitypara)
close(10)

if(iargc().eq.1)then
        CALL getarg(1,arg)
        READ(arg,*)argaline
        print *,'read in argaline = ',argaline
        else 
        argaline = 3.d0
endif


CALL INIT_STELLARDISK(n,domain)
CALL FindSpiral
dx = domain/dble(n)
dy = domain/dble(n)

ALLOCATE(density(2*n,2*n))
ALLOCATE(xcoord(2*n))
ALLOCATE(ycoord(2*n))

DO i = 1, n*2
        xcoord(i) =  dx*0.5d0 - domain + dble(i)*dx
        ycoord(i) =  dy*0.5d0 - domain + dble(i)*dx
ENDDO

ALLOCATE(shared_spiral.u(3,4*n))
ALLOCATE(shared_spiral.h1(4*n))
ALLOCATE(shared_spiral.phi1r(2*n))
ALLOCATE(shared_spiral.r(4*n))

shared_spiral = spiral
!$OMP PARALLEL SHARED(density,shared_spiral) PRIVATE(j,r,th,pi,pf,d)
spiral = shared_spiral
!$OMP DO PRIVATE(spiral)
DO i = 1, n*2
DO j = 1, n*2
        pf = (/xcoord(i),ycoord(j)/)
        pi = pf
        if(toproject)CALL projection(pi,pf)
        r  = sqrt(pi(1)**2+pi(2)**2)
        th = atan2(pi(2),pi(1))
        d  = sigma1(r,th)
        !===================
        !ignore d too high
        if(r.gt.12.d0)d = 0.d0
        !ignore inside r=1.26  
!       if(r.lt.2.26)d = d*(1.d0 - cos(r/2.26d0*pi_n/2.d0))
        !ignore d that is not exist during coordinate transformation
        if(isnan(d))d = 0.d0
        !ignore value below detection limit
!       if(abs(d).lt.limit)d = 0.d0
        density(i,j) = d
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 


!
!!!Find 2d Potential
 ALLOCATE(potential(2*n,2*n))
!DO i = 1, n*2
!DO j = 1, n*2
!        r = sqrt(xcoord(i)**2+ycoord(j)**2)
!        th = atan2(ycoord(j),xcoord(i))
!        potential(i,j) = phi1(r,th)
!ENDDO
!ENDDO

!!Find Force
ALLOCATE(force(n/4,n/4,2))
!DO i = 1,n/4
!DO j = 1, n/4
!        r = sqrt(xcoord(i*8)**2+ycoord(j*8)**2)
!        th = atan2(ycoord(j*8),xcoord(i*8))
!        call FindForce(force(i,j,:),r,th)
!ENDDO
!ENDDO

points(1,:) = (/0.0,10.636/)
points(2,:) = (/0.0,-10.636/)
points(3,:) = (/0.0,4.727/)
points(4,:) = (/0.0,-4.727/)
CALL dprojection(points)
CALL plotdensity(density,potential,force,n,domain)

DEALLOCATE(potential)
DEALLOCATE(xcoord)
DEALLOCATE(ycoord)
DEALLOCATE(density)
CALL ENDSTELLARDISK
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
