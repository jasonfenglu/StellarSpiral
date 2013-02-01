MODULE STELLARDISK
IMPLICIT NONE
DOUBLE PRECISION,PARAMETER::GravConst = 4.3d-6 
DOUBLE PRECISION,PARAMETER::g = 4.3d0
DOUBLE PRECISION,PARAMETER::pi=4.d0*atan(1.d0)
DOUBLE PRECISION          ::domain
!DOUBLE PRECISION          ::wr=57.102d0
!DOUBLE PRECISION          ::wi=-2.182d0
DOUBLE PRECISION,POINTER  ::para(:)=>null()
!Pspd without curF
DOUBLE PRECISION          ::wr=63.444d0
DOUBLE PRECISION          ::wi=-1.048d0
!solved wave
DOUBLE COMPLEX,ALLOCATABLE              ::u(:,:),h1(:),phi1r(:)
CONTAINS

SUBROUTINE INIT_STELLARDISK(n,ddomain)
IMPLICIT NONE
DOUBLE PRECISION        ::ddomain
INTEGER                 ::n
DOUBLE PRECISION,TARGET         ::stdpara(10)
!Halo
DOUBLE PRECISION                ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION                ::rb,Mb,gBulge,VBulge
!disk
DOUBLE PRECISION                ::dM,da,db,VDisk
!Toomre Q
DOUBLE PRECISION                ::Q,Qod,rq

Lh   = 2.8d0
rhoh = 4.0e7
Mb   = 10.0d7
rb   = 2.0d0
dM   = 7.0d10
da   = 2.7
db   = 0.3
Qod  = 1.d0
q    = 1.2d0
rq   = 2.8d0

stdpara = (/Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq/)
para=>stdpara

!!Allocate
ALLOCATE(u(3,4*n))
ALLOCATE(h1(4*n))
ALLOCATE(phi1r(2*n))

domain = ddomain

ENDSUBROUTINE INIT_STELLARDISK

SUBROUTINE CALSPIRAL(i,j)
IMPLICIT NONE
include 'omp_lib.h'
integer i,j

!find EigenFunction
CALL findu(u)
!find h1
CALL findh1(u,h1)
!!!Find phi1 along r
!CALL FindPhi1()
!CALL k3sqrtlog

ENDSUBROUTINE

SUBROUTINE k3sqrtlog
IMPLICIT NONE                           
INTEGER                                 ::i
DOUBLE PRECISION                        ::r
open(10,file='r-dep.dat')
DO i = 2, size(u,2),2
         r = real(u(1,i))
        write(10,'(5(1XE15.6))')real(u(1,i)),real(u(2,i)),real(h1(i))/snsd(r)**2*sigma0(r),real(phi1r(i/2)),real(h1(i))
enddo
close(10)
ENDSUBROUTINE

function ToomreQ(r)
DOUBLE PRECISION  Q,r,Qod,ToomreQ,rq

!Qod = 1.d0
!q   = 1.2d0
!rq  = 2.8d0
Qod = para(8)
q   = para(9)
rq  = para(10)
ToomreQ = Qod*(1.d0 + q*dexp(-r**2/rq**2))
endfunction

function nu(r)
IMPLICIT NONE
DOUBLE COMPLEX                  ::nu   
DOUBLE PRECISION                ::r
DOUBLE PRECISION                ::m 

m = 2.d0

nu = (dcmplx(wr,wi)-m*Omega(r))/kappa(r)

endfunction

function k3sqrt(r)
IMPLICIT NONE
DOUBLE COMPLEX            k3sqrt
DOUBLE PRECISION,INTENT(in)::r
DOUBLE PRECISION          ::rr


k3sqrt = (dcmplx(kappa(r)/snsd(r)))**2*(dcmplx(ToomreQ(r))**-2  &
         - 1.d0 + nu(r)**2)

!k3sqrt = (dcmplx(kappa(r)/snsd(r)))**2*(dcmplx(ToomreQ(r))**-2 &
!         - 1.d0 + nu(r)**2 + 0.25d0*curF(r)**2*ToomreQ(r)**2)
!print *,r,reappal(k3sqrt)
!print *,r,curf(r)

endfunction

function snsd(r)
IMPLICIT NONE
DOUBLE PRECISION          ::r,snsd

snsd = ToomreQ(r)*pi*g*sigma0(r)/kappa(r)

ENDFUNCTION

recursive function Sigma0(r)
IMPLICIT NONE
DOUBLE PRECISION          ::Sigma0,r
DOUBLE PRECISION          ::m,a,b
DOUBLE PRECISION          ::BOUND,EPSREL,EPSABS
DOUBLE PRECISION          ::ans
DOUBLE PRECISION          ::ABSERR
INTEGER                   ::NEVAL,IERR,LIMIT,LENW,LAST,INF
DOUBLE PRECISION,ALLOCATABLE ::WORK(:)
INTEGER,ALLOCATABLE       ::IWORK(:)

!M     = 7.0d10
!a     = 2.7
!b     = 0.3

M     = para(5)
a     = para(6)
b     = para(7)

BOUND = 0.d0
INF   = 2
EPSREL = 10d-10
EPSABS = 10d-10
LIMIT  = 100
LENW   = LIMIT*4+2
ALLOCATE(IWORK(LENW))
ALLOCATE(WORK(LENW))
CALL DQAGI(FUN,BOUND,INF,EPSABS,EPSREL,ANS,ABSERR,NEVAL,IERR,  &
           LIMIT,LENW,LAST,IWORK,WORK)
DEALLOCATE(WORK)
DEALLOCATE(IWORK)
Sigma0 = ans/10.d5
contains 
function FUN(z)
IMPLICIT NONE
DOUBLE PRECISION          ::fun,z
fun =  &
      (b**2*M/4.d0/pi)* &
      (a*r**2+(a+3.d0*sqrt(z**2+b**2))*(a+sqrt(z**2+b**2))**2)/ &
      (r**2+(a+sqrt(z**2+b**2))**2)**2.5/(z**2+b**2)**1.5

ENDFUNCTION

ENDFUNCTION

function kappa(r)
IMPLICIT NONE
DOUBLE PRECISION  kappa,r
DOUBLE PRECISION  dr
DOUBLE PRECISION  dOmega
kappa = sqrt(4.d0*Omega(r)**2*(1.d0+r/(2.d0*Omega(r))*dfunc(Omega,r)))
ENDFUNCTION

FUNCTION Omega(r)
IMPLICIT NONE
DOUBLE PRECISION          ::Omega
DOUBLE PRECISION,INTENT(IN)::r
!Halo
DOUBLE PRECISION          ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION          ::rb,Mb,gBulge,VBulge
!disk
DOUBLE PRECISION          ::dM,da,db,VDisk
!Halo
Lh   = para(1)
rhoh = para(2)
gHalo = 4.d0*Lh**2.d0*pi*rhoh*(r - Lh*atan(r/Lh))
gHalo = GravConst/(r**2)*gHalo
VHalo = sqrt(r*gHalo)

!Bulge
Mb   = para(3)
rb   = para(4)
gBulge = 4*pi*rb**3*Mb*(-r/sqrt(1+r**2/rb**2)/rb+asinh(r/rb))
gBulge = gBulge*GravConst/r**2
VBulge = sqrt(r*gBulge)

!Disk
dM     = para(5)
da     = para(6)
db     = para(7)
!VDisk  = sqrt(dfunc(pDisk,r)*r)
VDisk  = sqrt(dpDisk(r)*r)

Omega  = sqrt(VHalo**2+VBulge**2+VDisk**2)/r
CONTAINS
FUNCTION pDisk(r)
IMPLICIT NONE
DOUBLE PRECISION        ::pDisk,r
pDisk  = -GravConst*dM
pDisk  = pDisk/sqrt(r**2+(da+db)**2)
ENDFUNCTION

FUNCTION dpDisk(r)
IMPLICIT NONE
DOUBLE PRECISION        ::dpDisk
DOUBLE PRECISION        ::r
dpDisk = GravConst*dM*r
dpDisk = dpDisk/((da+db)**2+r**2)**1.5
ENDFUNCTION

ENDFUNCTION

function dfunc(func,r)
!
! Forward differential
!
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE PRECISION,PARAMETER      ::dr = 1.d-4,coe(3)=(/-1.5d0,2.d0,-0.5d0/)
DOUBLE PRECISION                ::dfunc,ans,funcs(3)
INTEGER                         ::i
interface 
        function func(x)
        DOUBLE PRECISION        ::func
        DOUBLE PRECISION        ::x
        ENDFUNCTION func
endinterface        

funcs(1) = func(r)
funcs(2) = func(r+dr)
funcs(3) = func(r+2.d0*dr)
ans = dot_product(funcs,coe)/dr
dfunc = ans

!ans = 0.d0
!ans =        -3.d0/2.d0*func(r)
!ans = ans +        2.d0*func(r+dr)
!ans = ans +  -1.d0/2.d0*func(r+2.d0*dr)
!dfunc = ans/dr
endfunction

function curf(r)
IMPLICIT NONE
DOUBLE PRECISION                ::curf
DOUBLE PRECISION                ::s,r,tmp
INTEGER                         ::m=2

s    = -r/Omega(r)*dfunc(Omega,r)
curf = 2*dble(m)*(pi*g*sigma0(r))/kappa(r)**2/r
curf = curf/sqrt(1.d0/s-0.5d0)

ENDFUNCTION

SUBROUTINE findu(u)
USE RK,only:rk4
IMPLICIT NONE
DOUBLE COMPLEX                  ::u(:,:)
DOUBLE COMPLEX                  ::ui(3)
DOUBLE PRECISION                ::a,b
!DOUBLE PRECISION                ::domain
INTEGER                         ::n

u = (1.d0,0.d0)*0.d0
n = size(u,2)
a = 0.001d0
b = 1.d0*domain
ui = (/a,1.d0,0.d0/)
CALL rk4(a,b,N,p,q,p,u,ui)
contains

RECURSIVE FUNCTION p(r)
IMPLICIT NONE
DOUBLE COMPLEX                  ::p
DOUBLE PRECISION,INTENT(IN)     ::r
p = (0.d0,0.d0)
ENDFUNCTION

RECURSIVE FUNCTION q(r)
IMPLICIT NONE
DOUBLE COMPLEX                  ::q
DOUBLE PRECISION,INTENT(IN)     ::r
q = k3sqrt(r)
!q = dsin(r)
ENDFUNCTION

ENDSUBROUTINE

SUBROUTINE findh1(u,h1)
IMPLICIT NONE
include 'omp_lib.h'
DOUBLE COMPLEX                  ::u(:,:),h1(:)
DOUBLE COMPLEX                  ::rad
DOUBLE PRECISION                ::h,r
INTEGER                         ::i,n
n = size(u,2)
do i = 1, N 
!find h1 by the interploted u
        r   = u(1,i)
        rad = sqrt(kappa(r)**2*(1.d0-nu(r)**2)/sigma0(r)/r)
        h1(i) = u(2,i)*rad*exp(-0.5d0*(0.d0,1.d0)*ExpPart(r))
enddo

CALL refineh1(u,h1)

CONTAINS
Function ExpPart(r)
        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(IN)     ::r
        DOUBLE PRECISION                ::rr
        DOUBLE PRECISION                ::ExpPart
        DOUBLE PRECISION                ::err = 10.d-8
        INTEGER                         ::IERR
        rr = r
        CALL DGAUS8(Sigma,0.d0,rr,ERR,ExpPart,IERR)

ENDFUNCTION ExpPart

ENDSUBROUTINE findh1

SUBROUTINE refineh1(u,h1)
IMPLICIT NONE
DOUBLE COMPLEX                  ::u(:,:),h1(:)
INTEGER                         ::i,f
DO i =1, size(u,2)
        if(abs(u(1,i)).gt.9.5d0)then
                f = i
                exit
        endif
enddo

DO i = f+1, size(u,2)
        h1(i) = h1(i)*exp((-abs(u(1,i))+9.5d0)/0.5d0)
enddo

ENDSUBROUTINE refineh1

SUBROUTINE FindPhi1()
IMPLICIT NONE
DOUBLE COMPLEX                  ::k(4)
!DOUBLE COMPLEX                  ::phi1r(:)
DOUBLE PRECISION                ::r,h
INTEGER                         ::i,j,l,n

n = size(u,2)
!!Solve ODE of phi from density by RK4
h = u(1,3)-u(1,1)
phi1r(1) = (0.d0,0.d0)
DO i = 2,N-2,2
        r = u(1,i-1)
        k(1) = h*dsimplifiedPoisson(r           ,phi1r(i/2)                ,h1(i-1))
        k(2) = h*dsimplifiedPoisson(r+h/2.d0    ,phi1r(i/2)+k(1)/2.d0      ,h1(i))
        k(3) = h*dsimplifiedPoisson(r+h/2.d0    ,phi1r(i/2)+k(2)/2.d0      ,h1(i))
        k(4) = h*dsimplifiedPoisson(r+h         ,phi1r(i/2)+k(3)           ,h1(i+1))
        phi1r(i/2+1) = phi1r(i/2) + (k(1)+2.d0*k(2)+2.d0*k(3)+k(4))/6.d0
ENDDO


ENDSUBROUTINE

FUNCTION dsimplifiedPoisson(r,phi,h)
IMPLICIT NONE
DOUBLE COMPLEX                  ::dsimplifiedPoisson
DOUBLE COMPLEX,INTENT(IN)       ::phi,h
DOUBLE PRECISION                ::r

dsimplifiedPoisson = -phi/(2.d0*r)+(0.d0,1.d0)*cmplx(Sigma(r),0)*h
!dsimplifiedPoisson = dsimplifiedPoisson + 3.75d0/(0.d0,1.d0)/Sigma(r)/r**2*phi

!!test case 
!!dsimplifiedPoisson =  -phi*r + (0.d0,1.d0)*r
!!bnd condition is phi = 1+i at r=0
!!solution is
!!phi = exp(-r**2/2)+ i
!!dsimplifiedPoisson = -phi*r
ENDFUNCTION

recursive function Sigma(r) RESULT(ans)
!This is NOT related to density
IMPLICIT NONE
DOUBLE PRECISION                ::ans
DOUBLE PRECISION,INTENT(IN)     ::r
ans = 2.d0*pi*G*sigma0(r)/snsd(r)**2
ENDFUNCTION

SUBROUTINE FindForce(force,r,th)
IMPLICIT NONE
DOUBLE PRECISION                ::force(2),r,th
DOUBLE COMPLEX                  ::hh1,phir
DOUBLE PRECISION                ::runit(2),thunit(2)
INTEGER                         ::i,n
!interploting u at non-grid point r
n = size(u,2)/2
do i = 1, n
        if(real(u(1,2*i)).gt.r)then
                exit
        endif
enddo
phir =(r-u(1,2*i-2))/(u(1,2*i)-u(1,2*i-2))*(phi1r(i)-phi1r(i-1)) + phi1r(i-1)
i = i*2
hh1 =(r-u(1,i-1))/(u(1,i)-u(1,i-1))*(h1(i)-h1(i-1)) + h1(i-1)

runit = (/cos(th),sin(th)/)
thunit = (/-sin(th),cos(th)/)

force = real(dsimplifiedPoisson(r,phir,hh1)*exp((0.d0,-2.d0)*th)*runit &
      + (0.d0,-2.d0)*phir/r*exp((0.d0,-2.d0)*th)*thunit)

ENDSUBROUTINE

FUNCTION sigma1(r,th)
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
DOUBLE PRECISION                ::sigma1
DOUBLE PRECISION,INTENT(IN)     ::r,th
DOUBLE COMPLEX                  ::uu,hh1
DOUBLE PRECISION                ::rad
INTEGER                         ::i,j,k,l,n
!interploting u at non-grid point r
n = size(u,2)
do i = 1, n
        if(real(u(1,i)).gt.r)then
                exit
        endif
enddo
!uu = (r-u(1,i-1))/(u(1,i)-u(1,i-1))*(u(2,i)-u(2,i-1)) + u(2,i-1)
hh1 =(r-u(1,i-1))/(u(1,i)-u(1,i-1))*(h1(i)-h1(i-1)) + h1(i-1)


!find density
sigma1 = real(hh1*sigma0(r)/snsd(r)**2*exp(-2.d0*th*(0.d0,1.d0)))


ENDFUNCTION

FUNCTION phi1(r,th)
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
DOUBLE PRECISION                ::phi1
DOUBLE PRECISION,INTENT(IN)     ::r,th
DOUBLE COMPLEX                  ::uu,hh1
DOUBLE PRECISION                ::rad
INTEGER                         ::i,j,k,l,n
!interploting u at non-grid point r
n = size(h1)
do i = 1, n*2
        if(real(u(1,2*i)).gt.r)then
                exit
        endif
enddo
!uu = (r-u(1,i-1))/(u(1,i)-u(1,i-1))*(u(2,i)-u(2,i-1)) + u(2,i-1)
hh1 =(r-u(1,2*i-2))/(u(1,2*i)-u(1,2*i-2))*(phi1r(i)-phi1r(i-1)) + phi1r(i-1)


!find potential
phi1 = real(hh1*exp(-2.d0*th*(0.d0,1.d0)))


ENDFUNCTION

SUBROUTINE ENDSTELLARDISK
DEALLOCATE(u)
DEALLOCATE(h1)
DEALLOCATE(phi1r)

ENDSUBROUTINE

function error()
IMPLICIT NONE          
DOUBLE COMPLEX          ::error
DOUBLE COMPLEX          ::uu(3)
DOUBLE PRECISION        ::h=10d-5   ,r
DOUBLE PRECISION        ::RE,AE,B,C,RR
INTEGER                 ::l,IFLAG

B = 6.d0
C = 12.d0
RR = 8.d0
RE = 1d-8
AE = 1d-8
call DFZERO(four21,B,C,RR,RE,AE,IFLAG)
r = B
do l = 1,size(u,2)
        if(real(u(1,l)).gt.r)then
                uu(:) = u(:,l)
                exit
        endif
enddo
error  = -(0.d0,1.d0)*sqrt(k3sqrt(r))
error  = error -                           &
0.5d0/sqrt(k3sqrt(r))*               &                
(sqrt(k3sqrt(r+h))-sqrt(k3sqrt(r-h)))/(2.d0*h)

error = uu(3)/uu(2) -error

CONTAINS
function four21(r)
IMPLICIT NONE
DOUBLE PRECISION                ::four21,r
four21 = (wr - 2.d0*Omega(r))/kappa(r) - 0.5d0
ENDFUNCTION
        
endfunction

SUBROUTINE single_grid(l,wri,wii)
!USE OMP_LIB
IMPLICIT NONE
include 'omp_lib.h'
type searchgrid_type
sequence
        DOUBLE PRECISION::coord(12,12,2)
        DOUBLE PRECISION::error(12,12)
endtype
type(searchgrid_type)           ::searchgrid,recvgrid
DOUBLE PRECISION                ::dr,wri,wii,di
INTEGER                         ::l,i,j,p(2)



dr = 1.d0/10.0d0**(l-1)
di = 0.5d0/10.0d0**(l-1)
!most left and upper grid
wri = wri +(-6.d0+0.5d0)*dr
wii = wii +(-6.d0+0.5d0)*dr

DO i = 1,12
        searchgrid%coord(:,i,2) = dble(i-1)*dr + wii
        searchgrid%coord(i,:,1) = dble(i-1)*dr + wri
enddo
DO j = 1,12
CALL INIT_STELLARDISK(100,20.d0)
DO i = 1,12
        wr = searchgrid%coord(i,j,1)
        wi = searchgrid%coord(i,j,2)
        CALL CALSPIRAL(i,j)
ENDDO
ENDDO
        
p = MINLOC(searchgrid%error(:,:))
i = p(1)
j = p(2)
wri = searchgrid%coord(i,j,1)
wii = searchgrid%coord(i,j,2)
!if(j.eq.1 .or. j.eq.12 .or. i.eq.1 .or. i.eq.12)l = l -1
!DO i = 1,12
!DO j = 1,12
!        print *,searchgrid%coord(i,j,:),searchgrid%error(i,j)
!ENDDO
!ENDDO

ENDSUBROUTINE

SUBROUTINE omp_single_grid(l,wri,wii,err)
!USE OMP_LIB
IMPLICIT NONE
include 'omp_lib.h'
type searchgrid_type
sequence
        DOUBLE PRECISION::coord(12,12,2)
        DOUBLE PRECISION::error(12,12)
endtype
type(searchgrid_type)           ::searchgrid,recvgrid
DOUBLE PRECISION                ::dr,wri,wii,di,err
INTEGER                         ::l,i,j,p(2)

INTEGER                 ::ipc


dr = 1.d0/10.0d0**(l-1)
di = 0.5d0/10.0d0**(l-1)
!most left and upper grid
wri = wri +(-6.d0+0.5d0)*dr
wii = wii +(-6.d0+0.5d0)*di
DO i = 1,12
!       searchgrid%coord(:,i,2) = dble(i-1)*dr + wii
!       searchgrid%coord(i,:,1) = dble(i-1)*dr + wri
        searchgrid%coord(:,i,2) =                wii
        searchgrid%coord(i,:,1) =                wri
enddo
CALL INIT_STELLARDISK(100,20.d0)
!$OMP PARALLEL DO PRIVATE(wi,wr,u,h1) SHARED(searchgrid)
!OMP DO        
DO i = 1,12
DO j = 1,12
        wr = searchgrid%coord(i,j,1)
        wi = searchgrid%coord(i,j,2)
!       print *,OMP_GET_THREAD_NUM(),wr,wi,j
        CALL CALSPIRAL(i,j)
!       searchgrid%error(i,j) = abs(error())
ENDDO
ENDDO
!OMP END DO
!$OMP END PARALLEL DO
CALL ENDSTELLARDISK
        
p = MINLOC(searchgrid%error(:,:))
i = p(1)
j = p(2)
wri = searchgrid%coord(i,j,1)
wii = searchgrid%coord(i,j,2)
err = searchgrid%error(i,j)
!if(j.eq.1 .or. j.eq.12 .or. i.eq.1 .or. i.eq.12)l = l -1
!DO i = 1,12
!DO j = 1,12
!        print *,searchgrid%coord(i,j,:),searchgrid%error(i,j)
!ENDDO
!ENDDO

ENDSUBROUTINE

SUBROUTINE findpspsd(wri,wii)
IMPLICIT NONE
DOUBLE PRECISION                ::wri,wii,err
INTEGER                         ::l
do l = 1,8
        CALL omp_single_grid(l,wri,wii,err)
!       print *,wri,wii,err
enddo

ENDSUBROUTINE 

SUBROUTINE TEMP
type    type_u
       DOUBLE COMPLEX,POINTER           ::u(:,:) 
       DOUBLE PRECISION,POINTER         ::rcoord(:)
       INTEGER                          ::length
endtype

type    type_h1
       DOUBLE COMPLEX,POINTER           ::h1(:) 
       DOUBLE PRECISION,POINTER         ::rcoord(:)
       INTEGER                          ::length
endtype

type    type_phi1
       DOUBLE COMPLEX,POINTER           ::h1(:) 
       DOUBLE PRECISION,POINTER         ::rcoord(:)
       INTEGER                          ::length
endtype

type    type_stellarspiral
        type(type_u)                    ::u
        type(type_h1)                   ::h1
        type(type_phi1)                 ::phi
        DOUBLE PRECISION                ::rmin,rmax
endtype     

type::Q_type
        DOUBLE PRECISION,POINTER        ::Qod
        DOUBLE PRECISION,POINTER        ::q
        DOUBLE PRECISION,POINTER        ::rq
endtype       
type::comp_type
        !Halo
        DOUBLE PRECISION,POINTER        ::Lh,rhoh,gHalo,VHalo
        !bulge
        DOUBLE PRECISION,POINTER        ::rb,Mb,gBulge,VBulge
        !disk
        DOUBLE PRECISION,POINTER        ::dM,da,db,VDisk
endtype
type(type_stellarspiral)                ::stellarspiral
!SUBROUTINE mpi_single_grid(l,wri,wii)
!IMPLICIT NONE
!include 'mpif.h'
!type searchgrid_type
!sequence
!        DOUBLE PRECISION::coord(12,12,2)
!        DOUBLE PRECISION::error(12,12)
!endtype
!type(searchgrid_type)           ::searchgrid,recvgrid
!DOUBLE PRECISION                ::dr,wri,wii,di
!INTEGER                         ::l,i,j,p(2)
!!!mpi
!INTEGER                         ::ierr,nprocs,myid,trunknum
!integer                         ::mpi_grid_type,dextent
!integer                         ::oldtypes(0:1), blockcounts(0:1), &
!                                  offsets(0:1),itag,istat
!dr = 1.d0/10.0d0**(l-1)
!di = 0.5d0/10.0d0**(l-1)
!!most left and upper grid
!wri = wri +(-6.d0+0.5d0)*dr
!wii = wii +(-6.d0+0.5d0)*dr
!
!DO i = 1,12
!        searchgrid%coord(:,i,2) = dble(i-1)*dr + wii
!        searchgrid%coord(i,:,1) = dble(i-1)*dr + wri
!enddo
!call MPI_INIT(ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
!call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
!
!call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,dextent,ierr)
!offsets(0) = 0
!oldtypes(0) = MPI_DOUBLE_PRECISION
!blockcounts(0) = 12*12*2
!offsets(1) = 12*12*2*dextent
!oldtypes(1) = MPI_DOUBLE_PRECISION 
!blockcounts(1) = 12*12
!call MPI_TYPE_STRUCT(2,blockcounts,offsets,oldtypes,mpi_grid_type,ierr)
!
!call MPI_TYPE_CONTIGUOUS(12*12*2,MPI_DOUBLE_PRECISION,mpi_grid_type,ierr)
!call MPI_TYPE_COMMIT(mpi_grid_type,ierr)
!
!trunknum = 12/nprocs + 1
!if(myid.ne.0)then
!!print *,myid,min0(trunknum*myid+trunknum,12)
!DO j = 1,12
!!DO i = 1,12
!DO i = 1+trunknum*myid,min0(trunknum*myid+trunknum,12)
!        wr = searchgrid%coord(i,j,1)
!        wi = searchgrid%coord(i,j,2)
!        !CALL INIT_STELLARDISK(10,20.d0)
!        !searchgrid%error(i,j) = abs(error())
!        !CALL ENDSTELLARDISK
!ENDDO
!ENDDO
!endif
!itag = 2000
!if(myid.ne.0)then
!        CALL MPI_SEND(searchgrid,1,mpi_grid_type,0,itag,MPI_COMM_WORLD,ierr)
!elseif(myid.eq.0)then
!        DO i = 1,nprocs-1
!        CALL MPI_RECV(recvgrid  ,1,mpi_grid_type,i,itag,MPI_COMM_WORLD,istat,ierr)
!        searchgrid%coord(1+trunknum*i:min0(trunknum*myid+trunknum,12),:,:) &
!       =  recvgrid%coord(1+trunknum*i:min0(trunknum*myid+trunknum,12),:,:) 
!        searchgrid%error(1+trunknum*i:min0(trunknum*myid+trunknum,12),:) &
!       =  recvgrid%error(1+trunknum*i:min0(trunknum*myid+trunknum,12),:) 
!        ENDDO
!endif
!CALL MPI_FINALIZE(ierr)
!        
!p = MINLOC(searchgrid%error(:,:))
!i = p(1)
!j = p(2)
!wri = searchgrid%coord(i,j,1)
!wii = searchgrid%coord(i,j,2)
!!if(j.eq.1 .or. j.eq.12 .or. i.eq.1 .or. i.eq.12)l = l -1
!DO i = 1,12
!DO j = 1,12
!        print *,searchgrid%coord(i,j,:),searchgrid%error(i,j)
!ENDDO
!ENDDO
!
!ENDSUBROUTINE
ENDSUBROUTINE

ENDMODULE STELLARDISK