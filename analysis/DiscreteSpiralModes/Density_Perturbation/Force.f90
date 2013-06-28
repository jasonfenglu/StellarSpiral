PROGRAM density1
USE PLOTTING
USE STELLARDISK_MODEL
USE &
STELLARDISK,only:FindSpiral,pi_n=>pi,sigma1,SpiralForce,GravConst,Omega,Sigma0,StellarOmega
USE projections,only:argaline
USE io
USE math
IMPLICIT NONE
INTEGER                         ::i,j,k,l
CHARACTER(len=32)               ::arg
CHARACTER(len=20),PARAMETER     ::hdfname='force/force.h5'
DOUBLE PRECISION                ::domain= 10.d0,dx,dy,r,th,pf(2),pi(2)
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::force(:,:,:)
DOUBLE PRECISION,ALLOCATABLE    ::fr(:,:),frsorted(:,:)
INTEGER,ALLOCATABLE             ::sortindex(:)
DOUBLE PRECISION                ::limit = 100.d0
DOUBLE PRECISION                ::d,co,fmax,amp
DOUBLE PRECISION                ::x,y
INTEGER,PARAMETER               ::n=512
type(typspiral)                 ::spiral
INTEGER                         ::ierr
namelist /forcenml/               amp 

!read in density plot related options
open(10,file='para.list')
read(10,nml=forcenml)
close(10)

CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)
dx = domain/dble(n)*2.d0
dy = domain/dble(n)*2.d0

ALLOCATE(density(n,n))
ALLOCATE(xcoord(n))
ALLOCATE(ycoord(n))
ALLOCATE(force(n,n,3))
ALLOCATE(fr(n**2,3))
ALLOCATE(frsorted(n**2,3))
ALLOCATE(sortindex(n**2))

DO i = 1, n
        xcoord(i) = 0.5d0*dx - domain + dble(i-1)*dx
        ycoord(i) = 0.5d0*dy - domain + dble(i-1)*dy
ENDDO

!!Filling density map
!$OMP PARALLEL SHARED(density,spiral) PRIVATE(j,r,th,pi,pf,d)
!$OMP DO 
DO i = 1, n
DO j = 1, n
        pf = (/xcoord(i),ycoord(j)/)
        pi = pf
        r  = sqrt(pi(1)**2+pi(2)**2)
        th = atan2(pi(2),pi(1))
        d  = sigma1(r,th,spiral)
!       d  = sigma0(r,spiral)
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

if(abs(spiral.error).gt.1d-5)then
        write(0,*)'!!!!!! wrong pspd:'
        write(0,*)'error:',abs(spiral.error)
        write(0,*)'pspd:',spiral.w
endif

Force = 0.d0
!!Filling Force Map
!$OMP PARALLEL SHARED(density,Force) PRIVATE(j,k,l,dx,dy,r)
!$OMP DO 
DO i = 1, n
DO j = 1, n
        DO k = 1, n
        DO l = 1, n
                IF((i.eq.k).and.(j.eq.l))cycle
                dx = xcoord(i) - xcoord(k)
                dy = ycoord(j) - xcoord(l)
                r = sqrt(dx**2+dy**2)
                Force(i,j,1) = Force(i,j,1) - density(k,l)*GravConst/r**3*dx
                Force(i,j,2) = Force(i,j,2) - density(k,l)*GravConst/r**3*dy
        ENDDO
        ENDDO
        Force(i,j,:) = Force(i,j,:)/dble(n**2)*domain**2*4.d0
        Force(i,j,3) = &
        sqrt(Force(i,j,1)**2+Force(i,j,2)**2)
        fr(i+(j-1)*n,1) = sqrt(xcoord(i)**2+ycoord(j)**2)
        fr(i+(j-1)*n,2) = Force(i,j,3)
        pf = (/xcoord(i),ycoord(j)/)
        pi = pf
        r  = sqrt(pi(1)**2+pi(2)**2)
        th = atan2(pi(2),pi(1))
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 


points(1,:) = (/0.0,10.636/)
points(2,:) = (/0.0,-10.636/)
points(3,:) = (/0.0,4.727/)
points(4,:) = (/0.0,-4.727/)


!!Sorting force in r-direction
CALL DPSORT(fr(:,1),n**2,sortindex,1,ierr)
DO i = 1, n**2
        j = sortindex(i)
        frsorted(i,:) = fr(j,:)
        r = frsorted(i,1)
        frsorted(i,3) = Omega(r,spiral)/r**2
ENDDO

!!Find Corrotation
CALL findco(co,spiral)
!!Find max force near co
DO i = 1, n**2
        if(abs(frsorted(i,1) - co)<1.d-3)then
                fmax = max(frsorted(i,2),fmax)
                r = frsorted(i,1)
        ENDIF
ENDDO
fmax = fmax*amp
200 FORMAT(6(G12.4,3X))
write(6,200)'max at','spiral f','centri f','percentage'
write(6,200)r,fmax,StellarOmega(r,spiral)**2*r,fmax/StellarOmega(r,spiral)**2/r*100.d0

!save fr and centrifugal force
CALL h5write(frsorted(:,1),n**2,hdfname,'r')
CALL h5write(frsorted(:,2),n**2,hdfname,'fr')
CALL h5write(frsorted(:,3),n**2,hdfname,'centri')

CALL h5write(xcoord,N,hdfname,'xcoord')
CALL h5write(ycoord,N,hdfname,'ycoord')
CALL h5write(force(:,:,1),N,N,hdfname,'fx')
CALL h5write(force(:,:,2),N,N,hdfname,'fy')
CALL h5write(force(:,:,3),N,N,hdfname,'f')

DEALLOCATE(fr)
DEALLOCATE(Force)
DEALLOCATE(xcoord)
DEALLOCATE(ycoord)
DEALLOCATE(density)
!CALL PhaseIntegrate
!CALL ENDSTELLARDISK
STOP

END PROGRAM

SUBROUTINE findco(co,spiral)
USE STELLARDISK
type(typspiral)                        ::spiral
DOUBLE PRECISION                        ::co
DOUBLE PRECISION                        ::B,C,R,RE,AE
INTEGER                                 ::IFALG

B = 0.d0
C = 9.d0
R = 5.d0
RE = 1d-7
AE = 1d-7

CALL DFZERO(fcorrotation,B,C,R,RE,AE,IFALG)
co = B
CONTAINS
FUNCTION fcorrotation(r)
DOUBLE PRECISION                        ::fcorrotation,r
        fcorrotation = Omega(r,spiral) - real(spiral.w)/2.d0
ENDFUNCTION

ENDSUBROUTINE
