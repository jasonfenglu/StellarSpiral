MODULE NUM
DOUBLE PRECISION,PARAMETER                      ::zerolimit=1.d-6
ENDMODULE NUM

MODULE STELLARDISK
USE PLOTTING
IMPLICIT NONE
DOUBLE PRECISION,PARAMETER::GravConst   = 4.3d-6 
DOUBLE PRECISION,PARAMETER::g           = 4.3d0
DOUBLE PRECISION,PARAMETER::pi          = 4.d0*atan(1.d0)
LOGICAL,PARAMETER         ::withf       = .true.
DOUBLE PRECISION,POINTER,SAVE           ::para(:)=>null()
DOUBLE PRECISION          ::wr          
DOUBLE PRECISION          ::wi          
DOUBLE PRECISION,TARGET ,SAVE,ALLOCATABLE::stdpara(:)
type   spiral_type
       DOUBLE COMPLEX,ALLOCATABLE       ::u(:,:) 
       DOUBLE COMPLEX,ALLOCATABLE       ::h1(:)
       DOUBLE COMPLEX,ALLOCATABLE       ::phi1r(:)
       DOUBLE PRECISION,ALLOCATABLE     ::r(:)
       DOUBLE PRECISION                 ::rmax
       DOUBLE PRECISION                 ::rmin
       DOUBLE PRECISION                 ::fortoone
       INTEGER                          ::N
       LOGICAL                          ::ucaled    = .false.
       LOGICAL                          ::h1caled   = .false.
       LOGICAL                          ::phi1rcaled= .false.
endtype
type(spiral_type)         ::spiral
!$OMP THREADPRIVATE(spiral,wr,wi,stdpara)
CONTAINS

SUBROUTINE INIT_STELLARDISK(n,domain)
USE OMP_LIB
IMPLICIT NONE
DOUBLE PRECISION                ::domain
INTEGER                         ::n
!Halo
DOUBLE PRECISION                ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION                ::rb,Mb,gBulge,VBulge
!disk
DOUBLE PRECISION                ::dM,da,db,VDisk
!Toomre Q
DOUBLE PRECISION                ::Q,Qod,rq
!Lau Disk
DOUBLE PRECISION                ::a1,a2,M1,M2
!pspd from readin
DOUBLE PRECISION                ::w(4)
!NAME LIST
namelist /paralist/ Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2,w


!$OMP CRITICAL
open(10,file='para.list')
read(10,nml=paralist)
close(10)
!$OMP END CRITICAL


!ALLOCATE(stdpara(14))
if(.not.allocated(stdpara))ALLOCATE(stdpara(14))
stdpara = (/Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2/)
para => stdpara


!Allocate
if(.not.allocated(spiral.u))ALLOCATE(spiral.u(3,4*n))
if(.not.allocated(spiral.h1))ALLOCATE(spiral.h1(4*n))
if(.not.allocated(spiral.phi1r))ALLOCATE(spiral.phi1r(2*n))
if(.not.allocated(spiral.r))ALLOCATE(spiral.r(4*n))


spiral%rmin = 0.d0 
spiral%rmax = 2.d0*domain
spiral%N    = 4*n

!!choose pspd
if(withf)then
        wr = w(1)
        wi = w(2)
else
        wr = w(3)
        wi = w(4)
endif

ENDSUBROUTINE INIT_STELLARDISK

SUBROUTINE READINSTDPARA
IMPLICIT NONE
!Halo
DOUBLE PRECISION                ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION                ::rb,Mb,gBulge,VBulge
!disk
DOUBLE PRECISION                ::dM,da,db,VDisk
!Toomre Q
DOUBLE PRECISION                ::Q,Qod,rq
!Lau Disk
DOUBLE PRECISION                ::a1,a2,M1,M2
!pspd from readin
DOUBLE PRECISION                ::w(4)
!NAME LIST
namelist /paralist/ Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2,w

!$OMP CRITICAL
open(10,file='para.list')
read(10,nml=paralist)
close(10)
!$OMP END CRITICAL

!ALLOCATE(stdpara(14))
if(.not.allocated(stdpara))ALLOCATE(stdpara(14))
stdpara = (/Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2/)

ENDSUBROUTINE

SUBROUTINE FindSpiral
IMPLICIT NONE
!find EigenFunction
CALL findu
!find h1
CALL findh1
!Find phi1 along r
CALL FindPhi1
!Save calculation results
CALL k3sqrtlog
ENDSUBROUTINE

SUBROUTINE k3sqrtlog
IMPLICIT NONE                           
INTEGER                                 ::i
DOUBLE PRECISION                        ::r
DOUBLE COMPLEX,ALLOCATABLE              ::u(:,:),h1(:),phi1r(:)
ALLOCATE(u(3,spiral.N))
ALLOCATE(h1(spiral.N))
ALLOCATE(phi1r(spiral.N))

u       = spiral.u
h1      = spiral.h1
phi1r   = spiral.phi1r
open(10,file='r-dep.dat')
DO i = 2, spiral.n,2
        r = spiral.r(i)
        write(10,'(6(1XE15.6))')spiral.r(i),abs(spiral.u(2,i)),abs(spiral.h1(i))/snsd(r)**2*sigma0(r),abs(spiral.phi1r(i/2)),abs(spiral.h1(i))
        !r, u, sigma1,potential1,h1
enddo
close(10)

CALL rlog

ENDSUBROUTINE

SUBROUTINE rlog
USE plotting,only:plotlog
IMPLICIT NONE
DOUBLE PRECISION,ALLOCATABLE            ::dat(:,:)
DOUBLE PRECISION                        ::r
INTEGER                                 ::i,m

!r, rotation, k3, Q
m = 4

ALLOCATE(dat(4,spiral.n))
dat(1,:) = spiral.r
do i = 1, spiral.n
        r=spiral.r(i)
        dat(2,i) = Omega(r)
        dat(3,i) = dble(k3sqrt(r))
        dat(4,i) = ToomreQ(r)
enddo

CALL plotlog(dat,m,spiral.n)
DEALLOCATE(dat)
ENDSUBROUTINE

function ToomreQ(r)
DOUBLE PRECISION  Q,r,Qod,ToomreQ,rq
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


if(withf)then
        k3sqrt = (dcmplx(kappa(r)/snsd(r)))**2*(dcmplx(ToomreQ(r))**-2 &
                 - 1.d0 + nu(r)**2 + 0.25d0*curF(r)**2*ToomreQ(r)**2)
else                 
        k3sqrt = (dcmplx(kappa(r)/snsd(r)))**2*(dcmplx(ToomreQ(r))**-2  &
                 - 1.d0 + nu(r)**2)
endif
!if(isnan(real(k3sqrt)))CALL XERMSG('k3sqrt','k3sqrt','k3sqrt is nan.',-98,2)
k3sqrt = k3sqrt
CALL CheckResult
CONTAINS

SUBROUTINE CheckResult
IMPLICIT NONE
CHARACTER(len=72)                       ::errormsg
write(errormsg,"(E10.3)")r
errormsg = trim(errormsg)
errormsg = 'k3sqrt nan.@r = '//errormsg
if(isnan(real(k3sqrt)))CALL XERMSG('k3sqrt','k3sqrt',errormsg,-98,2)
ENDSUBROUTINE

endfunction

function snsd(r)
IMPLICIT NONE
DOUBLE PRECISION          ::r,snsd
snsd = ToomreQ(r)*pi*GravConst*sigma0(r)/kappa(r)
ENDFUNCTION
 
function Sigma0(rr)
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)::rr
DOUBLE PRECISION          ::Sigma0,r
DOUBLE PRECISION          ::m,a,b
DOUBLE PRECISION          ::BOUND,EPSREL,EPSABS
DOUBLE PRECISION          ::ans
DOUBLE PRECISION          ::ABSERR
INTEGER                   ::NEVAL,IERR,LIMIT,LENW,LAST,INF
DOUBLE PRECISION,ALLOCATABLE ::WORK(:)
INTEGER,ALLOCATABLE       ::IWORK(:)

r = rr
!not using integral now
goto 10
M     = para(5)
a     = para(6)
b     = para(7)

BOUND = 0.d0
INF   = 2
EPSREL = 10d-12
EPSABS = 10d-15
LIMIT  = 2000
LENW   = LIMIT*4+2
ALLOCATE(WORK(LENW))
ALLOCATE(IWORK(LIMIT))
CALL DQAGI(FUN,BOUND,INF,EPSABS,EPSREL,ANS,ABSERR,NEVAL,IERR,  &
           LIMIT,LENW,LAST,IWORK,WORK)
DEALLOCATE(WORK)
DEALLOCATE(IWORK)
Sigma0 = ans/10.d5

10 sigma0 = LauDiskSigma0(r)

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

FUNCTION LauDiskSigma0(r)
IMPLICIT NONE
DOUBLE PRECISION                ::LauDiskSigma0
DOUBLE PRECISION                ::r
!Lau Disk
DOUBLE PRECISION                ::x1,x2
DOUBLE PRECISION                ::a1,a2,M1,M2
a1 = para(11)
a2 = para(12)
M1 = para(13)
M2 = para(14)
x1 = sqrt(1.d0+(r/a1)**2)
x2 = sqrt(1.d0+(r/a2)**2)

LauDiskSigma0 = 4.5d0*M1/pi/a1**2/x1**11
LauDiskSigma0 = LauDiskSigma0 - 4.5d0*M2/pi/a2**2/x2**11 

CONTAINS

FUNCTION H(x)
IMPLICIT NONE
DOUBLE PRECISION                ::H,x

H = 0.d0
H = 59.0625/x**11 + H
H = 26.25  /x**9  + H
H = 16.875 /x**7  + H
H = 11.25  /x**5  + H
H = 6.5625 /x**3  + H

ENDFUNCTION
ENDFUNCTION

function kappa(r)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION  kappa
DOUBLE PRECISION,INTENT(IN)             ::r
DOUBLE PRECISION  dr
DOUBLE PRECISION  dOmega
!if(r.lt.zerolimit)then
!        kappa = 2.d0*Omega(r)
!        return
!endif
kappa = sqrt(4.d0*Omega(r)**2*(1.d0+r/(2.d0*Omega(r))*dfunc(Omega,r)))
if(isnan(kappa))then
        print *,'kappa exception catched'
        print *,'r,Omega,dfunc(Omega,r):'
        print *,r,Omega(r),dfunc(Omega,r)
        CALL XERMSG('k3sqrt','kappa','kappa is nan.',-97,2)
endif       
ENDFUNCTION

FUNCTION Omega(r)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION          ::Omega
DOUBLE PRECISION,INTENT(IN)::r
!Halo
DOUBLE PRECISION          ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION          ::rb,Mb,gBulge,VBulge
!disk
DOUBLE PRECISION          ::dM,da,db
DOUBLE COMPLEX            ::VDisk
DOUBLE PRECISION          ::a1,a2,M1,M2

if(.not.associated(para))then
        print *,'para not init'
        stop
endif        


!Halo
Lh   = para(1)
rhoh = para(2)
!Bulge
Mb   = para(3)
rb   = para(4)
!Disk
a1 = para(11)
a2 = para(12)
M1 = para(13)
M2 = para(14)


!!Limit case when r->0:
if(r.lt.zerolimit)then
        Omega = 4.d0/3.d0*pi*GravConst*(MB+RHOH) &
              + 16.d0/105*GravConst*(M1/a1**3-M2/a2**3)*120.d0
        Omega = Omega**0.5
        return
endif


!Halo
Lh   = para(1)
rhoh = para(2)
gHalo = 4.d0*Lh**2.d0*pi*rhoh*(r - Lh*atan(r/Lh))/r**2
gHalo = GravConst*gHalo
VHalo = sqrt(r*gHalo)
!Bulge
Mb   = para(3)
rb   = para(4)
gBulge = (-r/sqrt(1.d0+r**2/rb**2)/rb+dasinh(r/rb))/r**2

gBulge = 4.d0*pi*rb**3*Mb*gBulge*GravConst
VBulge = sqrt(r*gBulge)

!!Disk
!dM     = para(5)
!da     = para(6)
!db     = para(7)
!!VDisk  = sqrt(dfunc(pDisk,r)*r)
!VDisk  = sqrt(dpDisk(r)*r)

!LauDisk 
VDisk  = VLauDisk(r)
Omega = sqrt(VHalo**2+VBulge**2+dble(VDisk**2))/r
!print *,r,V
!CALL CheckResult
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

SUBROUTINE CheckResult
IMPLICIT NONE
CHARACTER(len=72)                       ::errormsg
write(errormsg,"(E10.3)")r
errormsg = trim(errormsg)
errormsg = 'Omega nan.@r = '//errormsg
if(isnan(Omega))CALL XERMSG('k3sqrt','Omega',errormsg,-99,2)
ENDSUBROUTINE

FUNCTION VLauDisk(r)
IMPLICIT NONE
DOUBLE COMPLEX                  ::VLauDisk
DOUBLE PRECISION                ::r
DOUBLE PRECISION                ::x1,x2
DOUBLE PRECISION                ::a1,a2,M1,M2

a1 = para(11)
a2 = para(12)
M1 = para(13)
M2 = para(14)
x1 = sqrt(1.d0+(r/a1)**2)
x2 = sqrt(1.d0+(r/a2)**2)

VLauDisk = (M1*GravConst/a1**3*H(x1))
VLauDisk = VLauDisk - (M2*GravConst/a2**3*H(x2))
VLauDisk = sqrt(16.d0/105.d0*VLauDisk)
VLauDisk = VLauDisk*r
ENDFUNCTION

FUNCTION H(x)
IMPLICIT NONE
DOUBLE PRECISION                ::H,x

H = 0.d0
H = 59.0625/x**11 + H
H = 26.25  /x**9  + H
H = 16.875 /x**7  + H
H = 11.25  /x**5  + H
H = 6.5625 /x**3  + H

ENDFUNCTION

ENDFUNCTION

function dfunc(func,r)
!
! Forward differential
!
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE PRECISION                ::dr
DOUBLE PRECISION,PARAMETER      ::coe(3)=(/-1.5d0,2.d0,-0.5d0/)
DOUBLE PRECISION                ::dfunc,ans,funcs(3)
INTEGER                         ::i
interface 
        function func(x)
        DOUBLE PRECISION        ::func
        DOUBLE PRECISION        ::x
        ENDFUNCTION func
endinterface        

dr = epsilon(r)**0.3*max(r,epsilon(0d0))
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
if(isnan(dfunc))CALL XERMSG('k3sqrt','dfunc','dfunc is nan.',-94,0)
endfunction
 
function curf(r)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION                ::curf
DOUBLE PRECISION                ::r,tmp
INTEGER                         ::m=2

if(r.lt.zerolimit)then
        curf = &
        2.d0*dble(m)*(pi*GravConst*sigma0(0.d0))/kappa(0.d0)**2*sqrt(-2.d0*Omega2()/(Omega(0.d0)))
        return
elseif(r.lt.3.d-2)then
        curf = &
        2.d0*dble(m)*(pi*GravConst*sigma0(r))/kappa(r)**2 &
        *sqrt(-(2.d0*r*Omega2())/(r*Omega(r)))
else
        curf = &
        2.d0*dble(m)*(pi*GravConst*sigma0(r))/kappa(r)**2*sqrt(-dfunc(Omega,r)/(r*Omega(r)))
endif

CALL CheckResult
CONTAINS
FUNCTION Omega2
IMPLICIT NONE
DOUBLE PRECISION          ::Omega2
!Halo
DOUBLE PRECISION          ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION          ::rb,Mb,gBulge,VBulge
!disk
DOUBLE PRECISION          ::dM,da,db
DOUBLE COMPLEX            ::VDisk
DOUBLE PRECISION          ::a1,a2,M1,M2

if(.not.associated(para))then
        print *,'para not init'
        stop
endif        


!Halo
Lh   = para(1)
rhoh = para(2)
!Bulge
Mb   = para(3)
rb   = para(4)
!Disk
a1 = para(11)
a2 = para(12)
M1 = para(13)
M2 = para(14)

Omega2 = 4.d0*pi*GravConst*rhoh/5.d0/Lh**2 + 6.d0/5.d0*pi*GravConst*Mb/rb**2 &
       + 8.d0*GravConst/105.d0*(M1/a1**5-M2/a2**5)*1080.d0
Omega2 = -Omega2/2.d0/Omega(0.d0)
if(Omega2.gt.0.d0)then
        write(0,*)'Omega2 is larger than 0. curf will not exist. Omega2=',Omega2
        stop
endif

ENDFUNCTION

SUBROUTINE CheckResult
IMPLICIT NONE
CHARACTER(len=72)                       ::errormsg
errormsg = trim(errormsg)
errormsg = 'curf real nan.'
if(isnan(real(curf)))then
        print *,'curf exception'
        print *,r,Omega(r),dfunc(Omega,r)
        CALL XERMSG('k3sqrt','curf',errormsg,-96,2)
endif
ENDSUBROUTINE
ENDFUNCTION
 
SUBROUTINE findu
USE RK,only:rk4
IMPLICIT NONE
DOUBLE COMPLEX                  ::ui(3)
DOUBLE PRECISION                ::a,b
INTEGER                         ::n
LOGICAL                         ::ufzero
namelist /BND/ ufzero
spiral.u = (1.d0,0.d0)*0.d0
a = spiral.rmin
b = spiral.rmax
!$OMP CRITICAL
open(10,file='para.list')
read(10,nml=BND)
close(10)
!$OMP END CRITICAL
if(ufzero)then
        ui = (/dcmplx(a),dcmplx(0.d0),2.d0*sqrt(-q(0.d0))/)
else
        ui = (/dcmplx(a),dcmplx(1.d0,0.d0),dcmplx(0.d0)/)
endif
CALL rk4(a,b,spiral.n,p,q,p,spiral.u,ui)
spiral.ucaled = .true.
spiral.r = real(spiral.u(1,:))
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
!q =  (1)
!q = dcmplx(1.d0,r**2)
ENDFUNCTION

ENDSUBROUTINE

SUBROUTINE findh1()
IMPLICIT NONE
DOUBLE COMPLEX                  ::u
DOUBLE COMPLEX                  ::rad
DOUBLE PRECISION                ::h,r,expp
INTEGER                         ::i,n

do i = 1, spiral.n
        r   = spiral.r(i)
        u   = spiral.u(2,i)
        rad = sqrt(kappa(r)**2*(1.d0-nu(r)**2)/sigma0(r)/r)
        expp= ExpPart(r)
        spiral.h1(i) = u*rad*exp(-0.5d0*(0.d0,1.d0)*ExpPart(r))
!       spiral.h1(i) = u*rad*exp(-0.5d0*(0.d0,1.d0)*cmplx(expp,0.d0))
enddo

!CALL refineh1
spiral.h1caled = .true.

CONTAINS

Function ExpPart(r) result(ans)
        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(IN)     ::r
        DOUBLE PRECISION                ::a,rr
        DOUBLE PRECISION                ::ans
        DOUBLE PRECISION                ::ferr = 1.d-20
        INTEGER                         ::IERR,K=6000
        a = 1.d-6
        rr = r
        IERR = 0
        if(rr.eq.0.d0)then
                ans = 0.d0
        else
                CALL DGAUS8(Sigma,a,rr,fERR,ans,IERR)
!               CALL DQNC79(Sigma,a,rr,fERR,ans,IERR,K)
        endif
ENDFUNCTION ExpPart

SUBROUTINE refineh1
IMPLICIT NONE
DOUBLE PRECISION,PARAMETER      ::cutoff = 9.5
INTEGER                         ::i,f
DO i =1, spiral.n
        if(spiral.r(i).gt.cutoff)then
                f = i
                exit
        endif
enddo

DO i = f+1, spiral.n
        spiral.h1(i) = spiral.h1(i)*exp((-spiral.r(i)+cutoff)/1.0d-1)
enddo

ENDSUBROUTINE refineh1
ENDSUBROUTINE findh1

SUBROUTINE FindPhi1()
IMPLICIT NONE
DOUBLE COMPLEX                  ::k(4)
DOUBLE PRECISION                ::r,h
INTEGER                         ::i,j,l,n

n = spiral.n
!!Solve ODE of phi from density by RK4
h = spiral.u(1,3)-spiral.u(1,1)
spiral.phi1r(1) = (0.d0,0.d0)
DO i = 2,N-2,2
        r = spiral.r(i-1)
        k(1) = h*dsimplifiedPoisson(r           ,spiral.phi1r(i/2)                ,spiral.h1(i-1))
        k(2) = h*dsimplifiedPoisson(r+h/2.d0    ,spiral.phi1r(i/2)+k(1)/2.d0      ,spiral.h1(i))
        k(3) = h*dsimplifiedPoisson(r+h/2.d0    ,spiral.phi1r(i/2)+k(2)/2.d0      ,spiral.h1(i))
        k(4) = h*dsimplifiedPoisson(r+h         ,spiral.phi1r(i/2)+k(3)           ,spiral.h1(i+1))
        spiral.phi1r(i/2+1) = spiral.phi1r(i/2) + (k(1)+2.d0*k(2)+2.d0*k(3)+k(4))/6.d0
ENDDO


ENDSUBROUTINE

FUNCTION dsimplifiedPoisson(r,phi,h)
IMPLICIT NONE
DOUBLE COMPLEX                  ::dsimplifiedPoisson
DOUBLE COMPLEX,INTENT(IN)       ::phi,h
DOUBLE PRECISION                ::r

dsimplifiedPoisson = -phi/(2.d0*r)+(0.d0,1.d0)*cmplx(Sigma(r),0)*h
dsimplifiedPoisson = dsimplifiedPoisson + 3.75d0/(0.d0,1.d0)/Sigma(r)/r**2*phi

!!test case 
!!dsimplifiedPoisson =  -phi*r + (0.d0,1.d0)*r
!!bnd condition is phi = 1+i at r=0
!!solution is
!!phi = exp(-r**2/2)+ i
!!dsimplifiedPoisson = -phi*r
ENDFUNCTION

function Sigma(r) RESULT(ans)
!This is NOT related to density
IMPLICIT NONE
DOUBLE PRECISION                ::ans
DOUBLE PRECISION,INTENT(IN)     ::r
ans = 2.d0*pi*GravConst*sigma0(r)/snsd(r)**2
ENDFUNCTION

SUBROUTINE FindForce(force,r,th)
IMPLICIT NONE
DOUBLE PRECISION                ::force(2),r,th
DOUBLE COMPLEX                  ::hh1,phir
DOUBLE PRECISION                ::runit(2),thunit(2)
INTEGER                         ::i,n
!interploting u at non-grid point r
n = spiral.n/2
do i = 1, n
        if(spiral.r(2*i).gt.r)then
                exit
        endif
enddo
phir =(r-spiral.u(1,2*i-2))/(spiral.u(1,2*i)-spiral.u(1,2*i-2))*(spiral.phi1r(i)-spiral.phi1r(i-1)) + spiral.phi1r(i-1)
i = i*2
hh1 =(r-spiral.u(1,i-1))/(spiral.u(1,i)-spiral.u(1,i-1))*(spiral.h1(i)-spiral.h1(i-1)) + spiral.h1(i-1)

runit = (/cos(th),sin(th)/)
thunit = (/-sin(th),cos(th)/)

force = real(dsimplifiedPoisson(r,phir,hh1)*exp((0.d0,-2.d0)*th)*runit &
      + (0.d0,-2.d0)*phir/r*exp((0.d0,-2.d0)*th)*thunit)

ENDSUBROUTINE

FUNCTION sigma1(r,th)
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
DOUBLE COMPLEX                  ::uu,hh1
DOUBLE PRECISION                ::sigma1
DOUBLE PRECISION,INTENT(IN)     ::r,th
DOUBLE PRECISION                ::rad
INTEGER                         ::i,j,k,l,n
!interploting u at non-grid point r
n = spiral.n
do i = 1, n
        if(spiral.r(i).gt.r)then
                hh1 =(r-spiral.r(i-1))/(spiral.r(i)-spiral.r(i-1))*(spiral.h1(i)-spiral.h1(i-1)) + spiral.h1(i-1)
                exit
        endif
enddo
!uu = (r-u(1,i-1))/(u(1,i)-u(1,i-1))*(u(2,i)-u(2,i-1)) + u(2,i-1)


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
n = spiral.n
do i = 1, n*2
        if(spiral.r(2*i).gt.r)then
                exit
        endif
enddo
!uu = (r-u(1,i-1))/(u(1,i)-u(1,i-1))*(u(2,i)-u(2,i-1)) + u(2,i-1)
hh1 =(r-spiral.u(1,2*i-2))/(spiral.u(1,2*i)-spiral.u(1,2*i-2))*(spiral.phi1r(i)-spiral.phi1r(i-1)) + spiral.phi1r(i-1)


!find potential
phi1 = real(hh1*exp(-2.d0*th*(0.d0,1.d0)))


ENDFUNCTION

SUBROUTINE ENDSTELLARDISK
DEALLOCATE(spiral.u)
DEALLOCATE(spiral.h1)
DEALLOCATE(spiral.phi1r)
DEALLOCATE(spiral.r)
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
spiral.fortoone = r
do l = 1,spiral.n
        if(spiral.r(l).gt.r)then
                uu(:) = spiral.u(:,l)
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

SUBROUTINE TEMP
!type    type_u
!       DOUBLE COMPLEX,POINTER           ::u(:,:) 
!       DOUBLE PRECISION,POINTER         ::rcoord(:)
!       INTEGER                          ::length
!endtype
!
!type    type_h1
!       DOUBLE COMPLEX,POINTER           ::h1(:) 
!       DOUBLE PRECISION,POINTER         ::rcoord(:)
!       INTEGER                          ::length
!endtype
!
!type    type_phi1
!       DOUBLE COMPLEX,POINTER           ::h1(:) 
!       DOUBLE PRECISION,POINTER         ::rcoord(:)
!       INTEGER                          ::length
!endtype
!
!type    type_stellarspiral
!        type(type_u)                    ::u
!        type(type_h1)                   ::h1
!        type(type_phi1)                 ::phi
!        DOUBLE PRECISION                ::rmin,rmax
!endtype     
!
!type::Q_type
!        DOUBLE PRECISION,POINTER        ::Qod
!        DOUBLE PRECISION,POINTER        ::q
!        DOUBLE PRECISION,POINTER        ::rq
!endtype       
!type::comp_type
!        !Halo
!        DOUBLE PRECISION,POINTER        ::Lh,rhoh,gHalo,VHalo
!        !bulge
!        DOUBLE PRECISION,POINTER        ::rb,Mb,gBulge,VBulge
!        !disk
!        DOUBLE PRECISION,POINTER        ::dM,da,db,VDisk
!endtype
!type(type_stellarspiral)                ::stellarspiral
!!SUBROUTINE mpi_single_grid(l,wri,wii)
!!IMPLICIT NONE
!!include 'mpif.h'
!!type searchgrid_type
!!sequence
!!        DOUBLE PRECISION::coord(12,12,2)
!!        DOUBLE PRECISION::error(12,12)
!!endtype
!!type(searchgrid_type)           ::searchgrid,recvgrid
!!DOUBLE PRECISION                ::dr,wri,wii,di
!!INTEGER                         ::l,i,j,p(2)
!!!!mpi
!!INTEGER                         ::ierr,nprocs,myid,trunknum
!!integer                         ::mpi_grid_type,dextent
!!integer                         ::oldtypes(0:1), blockcounts(0:1), &
!!                                  offsets(0:1),itag,istat
!!dr = 1.d0/10.0d0**(l-1)
!!di = 0.5d0/10.0d0**(l-1)
!!!most left and upper grid
!!wri = wri +(-6.d0+0.5d0)*dr
!!wii = wii +(-6.d0+0.5d0)*dr
!!
!!DO i = 1,12
!!        searchgrid%coord(:,i,2) = dble(i-1)*dr + wii
!!        searchgrid%coord(i,:,1) = dble(i-1)*dr + wri
!!enddo
!!call MPI_INIT(ierr)
!!call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
!!call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
!!
!!call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,dextent,ierr)
!!offsets(0) = 0
!!oldtypes(0) = MPI_DOUBLE_PRECISION
!!blockcounts(0) = 12*12*2
!!offsets(1) = 12*12*2*dextent
!!oldtypes(1) = MPI_DOUBLE_PRECISION 
!!blockcounts(1) = 12*12
!!call MPI_TYPE_STRUCT(2,blockcounts,offsets,oldtypes,mpi_grid_type,ierr)
!!
!!call MPI_TYPE_CONTIGUOUS(12*12*2,MPI_DOUBLE_PRECISION,mpi_grid_type,ierr)
!!call MPI_TYPE_COMMIT(mpi_grid_type,ierr)
!!
!!trunknum = 12/nprocs + 1
!!if(myid.ne.0)then
!!!print *,myid,min0(trunknum*myid+trunknum,12)
!!DO j = 1,12
!!!DO i = 1,12
!!DO i = 1+trunknum*myid,min0(trunknum*myid+trunknum,12)
!!        wr = searchgrid%coord(i,j,1)
!!        wi = searchgrid%coord(i,j,2)
!!        !CALL INIT_STELLARDISK(10,20.d0)
!!        !searchgrid%error(i,j) = abs(error())
!!        !CALL ENDSTELLARDISK
!!ENDDO
!!ENDDO
!!endif
!!itag = 2000
!!if(myid.ne.0)then
!!        CALL MPI_SEND(searchgrid,1,mpi_grid_type,0,itag,MPI_COMM_WORLD,ierr)
!!elseif(myid.eq.0)then
!!        DO i = 1,nprocs-1
!!        CALL MPI_RECV(recvgrid  ,1,mpi_grid_type,i,itag,MPI_COMM_WORLD,istat,ierr)
!!        searchgrid%coord(1+trunknum*i:min0(trunknum*myid+trunknum,12),:,:) &
!!       =  recvgrid%coord(1+trunknum*i:min0(trunknum*myid+trunknum,12),:,:) 
!!        searchgrid%error(1+trunknum*i:min0(trunknum*myid+trunknum,12),:) &
!!       =  recvgrid%error(1+trunknum*i:min0(trunknum*myid+trunknum,12),:) 
!!        ENDDO
!!endif
!!CALL MPI_FINALIZE(ierr)
!!        
!!p = MINLOC(searchgrid%error(:,:))
!!i = p(1)
!!j = p(2)
!!wri = searchgrid%coord(i,j,1)
!!wii = searchgrid%coord(i,j,2)
!!!if(j.eq.1 .or. j.eq.12 .or. i.eq.1 .or. i.eq.12)l = l -1
!!DO i = 1,12
!!DO j = 1,12
!!        print *,searchgrid%coord(i,j,:),searchgrid%error(i,j)
!!ENDDO
!!ENDDO
!!
!!ENDSUBROUTINE
ENDSUBROUTINE

ENDMODULE STELLARDISK
