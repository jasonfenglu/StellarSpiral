PROGRAM density1
USE PLOTTING
USE STELLARDISK_MODEL
USE &
STELLARDISK,only:FindSpiral,pi_n=>pi,sigma1,phi1,FindPhi1,SpiralForce,GravConst,Omega,intplt,Sigma0,StellarOmega
USE projections,only:argaline
IMPLICIT NONE
INTEGER                         ::i,j,k,l
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::domain= 10.d0,dx,dy,r,th,pf(2),pi(2)
DOUBLE PRECISION,ALLOCATABLE    ::density(:,:),xcoord(:),ycoord(:)
DOUBLE PRECISION,ALLOCATABLE    ::force(:,:,:)
DOUBLE PRECISION,ALLOCATABLE    ::fr(:,:),frsorted(:,:)
INTEGER,ALLOCATABLE             ::sortindex(:)
DOUBLE PRECISION                ::limit = 100.d0
DOUBLE PRECISION                ::d,co,fmax,fratio
INTEGER,PARAMETER               ::n=200
type(typspiral)                 ::spiral
LOGICAL                         ::toproject
INTEGER                         ::ierr
namelist /densitypara/ toproject
namelist /forcenml/                fratio 

!read in density plot related options
open(10,file='para.list')
read(10,nml=densitypara)
read(10,nml=forcenml)
close(10)

CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)
CALL FindPhi1(spiral)
dx = domain/dble(n)
dy = domain/dble(n)

ALLOCATE(density(2*n,2*n))
ALLOCATE(xcoord(2*n))
ALLOCATE(ycoord(2*n))
ALLOCATE(force(2*n,2*n,3))
ALLOCATE(fr(4*n**2,3))
ALLOCATE(frsorted(4*n**2,3))
ALLOCATE(sortindex(4*n**2))

DO i = 1, n*2
        xcoord(i) =  dx*0.5d0 - domain + dble(i)*dx
        ycoord(i) =  dy*0.5d0 - domain + dble(i)*dx
ENDDO

!!Filling density map
!$OMP PARALLEL SHARED(density,spiral) PRIVATE(j,r,th,pi,pf,d)
!$OMP DO 
DO i = 1, n*2
DO j = 1, n*2
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
DO i = 1, n*2
DO j = 1, n*2
        !x-direction
        DO k = 1, n*2
        DO l = 1, n*2
                IF((i.ne.k).and.(j.ne.l))THEN
                        dx = xcoord(i) - xcoord(k)
                        dy = ycoord(j) - xcoord(l)
                        r = sqrt(dx**2+dy**2)
                        Force(i,j,:) = Force(i,j,:) - density(k,l)*GravConst/r**3*(/dx,dy/)
                ENDIF
        ENDDO
        ENDDO
        Force(i,j,3) = &
        sqrt(Force(i,j,1)**2+Force(i,j,2)**2)/dble(n**2)*domain**2
        fr(i+(j-1)*2*n,1) = sqrt(xcoord(i)**2+ycoord(j)**2)
        fr(i+(j-1)*2*n,2) = Force(i,j,3)
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 


points(1,:) = (/0.0,10.636/)
points(2,:) = (/0.0,-10.636/)
points(3,:) = (/0.0,4.727/)
points(4,:) = (/0.0,-4.727/)

CALL plotforce(Force(:,:,3),Force(:,:,2),n,domain)

!!Sorting force in r-direction
CALL DPSORT(fr(:,1),4*n**2,sortindex,1,ierr)
DO i = 1, 4*n**2
        j = sortindex(i)
        frsorted(i,:) = fr(j,:)
ENDDO

!!Find Corrotation
CALL findco(co,spiral)
!!Find max force near co
DO i = 1, 4*n**2
        if((abs(frsorted(i,1) - co)<1.d0).and.(frsorted(i,2)>fmax))then
                fmax = frsorted(i,2)
                r = frsorted(i,1)
        ENDIF
ENDDO
200 FORMAT(6(G12.4,3X))
write(6,200)'corrotation','spiral f','centri f','percentage','factor'
write(6,200)r,fmax,StellarOmega(r,spiral)**2*r,fmax/StellarOmega(r,spiral)**2/r*100.d0,fratio/fmax*StellarOmega(r,spiral)**2*r/100.d0


!!!Print fr and centirfugal force
open(10,file='force.log')
DO i = 1, 4*n**2
        r = frsorted(i,1)
        write(10,*)r,frsorted(i,2),Omega(r,spiral)/r**2
ENDDO
close(10)


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
