PROGRAM caldensity
USE PLOTTING
USE STELLARDISK
IMPLICIT NONE
INTEGER                         ::i,j,k
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::domain= 12.d0,dx,dy,r,th
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::potential(:,:)
DOUBLE PRECISION,ALLOCATABLE    ::force(:,:,:)
INTEGER,PARAMETER               ::n=500
type(spiral_type)               ::shared_spiral

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
!$OMP PARALLEL SHARED(density,shared_spiral) PRIVATE(j,r,th)
spiral = shared_spiral
!$OMP DO PRIVATE(spiral)
DO i = 1, n*2
DO j = 1, n*2
        r = sqrt(xcoord(i)**2+ycoord(j)**2)
        th = atan2(ycoord(j),xcoord(i))
        density(i,j) = sigma1(r,th)
!       density(i,j) = sigma0(r)
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 

open(10,file='r-dep.dat')
DO i = 2, spiral.n,2
        r = spiral.r(i)
        write(10,'(5(1XE15.6))')spiral.r(i),real(spiral.u(2,i)),real(spiral.h1(i))/snsd(r)**2*sigma0(r),real(spiral.phi1r(i/2)),real(spiral.h1(i))
        !r, u, sigma1,potential1
enddo
close(10)


!!Find 2d Potential
ALLOCATE(potential(2*n,2*n))
DO i = 1, n*2
DO j = 1, n*2
        r = sqrt(xcoord(i)**2+ycoord(j)**2)
        th = atan2(ycoord(j),xcoord(i))
        potential(i,j) = phi1(r,th)
ENDDO
ENDDO

!!Find Force
ALLOCATE(force(n/4,n/4,2))
!DO i = 1,n/4
!DO j = 1, n/4
!        r = sqrt(xcoord(i*8)**2+ycoord(j*8)**2)
!        th = atan2(ycoord(j*8),xcoord(i*8))
!        call FindForce(force(i,j,:),r,th)
!ENDDO
!ENDDO

CALL plotdensity(density,potential,force,n,domain)
DEALLOCATE(potential)
DEALLOCATE(xcoord)
DEALLOCATE(ycoord)
DEALLOCATE(density)
CALL ENDSTELLARDISK
STOP


END PROGRAM

