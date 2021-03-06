PROGRAM calforce
USE PLOTTING
USE STELLARDISK
IMPLICIT NONE
INTEGER                         ::i,j,k
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::domain= 12.d0,dx,dy,r,th
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::potential(:,:)
DOUBLE PRECISION,ALLOCATABLE    ::force(:,:,:)
DOUBLE PRECISION                ::tmp
INTEGER,PARAMETER               ::n=100


CALL INIT_STELLARDISK(n,domain)

dx = domain/dble(n)
dy = domain/dble(n)

ALLOCATE(density(2*n,2*n))
ALLOCATE(xcoord(2*n))
ALLOCATE(ycoord(2*n))

DO i = 1, n*2
        xcoord(i) =  dx*0.5d0 - domain + dble(i)*dx
        ycoord(i) =  dy*0.5d0 - domain + dble(i)*dx
ENDDO

!DO i = 1, n*2
!DO j = 1, n*2
!        r = sqrt(xcoord(i)**2+ycoord(j)**2)
!        th = atan2(ycoord(j),xcoord(i))
!        density(i,j) = sigma1(r,th)
!!       density(i,j) = sigma0(r)
!ENDDO
!ENDDO


open(10,file='r-dep.dat')
DO i = 2, n*4,2
         r = spiral.r(i)
        write(10,'(5(1XE15.6))')spiral.r(i),real(spiral.u(2,i)),real(spiral.h1(i))/snsd(r)**2*sigma0(r),real(spiral.phi1r(i/2)),real(spiral.h1(i))
enddo
close(10)

!!Find 2d Potential
ALLOCATE(potential(2*n,2*n))
!DO i = 1, n*2
!DO j = 1, n*2
!        r = sqrt(xcoord(i)**2+ycoord(j)**2)
!        th = atan2(ycoord(j),xcoord(i))
!        potential(i,j) = phi1(r,th)
!ENDDO
!ENDDO

!!Find Force
ALLOCATE(force(2*n,2*n,2))
DO i = 1, n*2
DO j = 1, n*2
        r = sqrt(xcoord(i)**2+ycoord(j)**2)
        th = atan2(ycoord(j),xcoord(i))
        call FindForce(force(i,j,:),r,th)
        tmp = sqrt(force(i,j,1)**2+force(i,j,2)**2)
        write(*,'(5(1XE15.6))')r,tmp,tmp/(Omega(r)**2*r)
ENDDO
ENDDO

!CALL plot2d(density,potential,force,n,domain)
DEALLOCATE(potential)
DEALLOCATE(xcoord)
DEALLOCATE(ycoord)
DEALLOCATE(density)
CALL ENDSTELLARDISK
STOP


END PROGRAM
