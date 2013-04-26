MODULE NUM
DOUBLE PRECISION,PARAMETER                      ::zerolimit=1.d-6
ENDMODULE NUM

MODULE STELLARDISK_MODEL
type   typgalaxy_para
       DOUBLE PRECISION                 ::para(14)
       LOGICAL                          ::inited=.false.
       CONTAINS 
       PROCEDURE::printpara             =>para_print
       PROCEDURE,NOPASS::readstd        =>READINSTDPARA
       PROCEDURE,NOPASS::cpfromstd      =>cpfromstd
endtype
type(typgalaxy_para),TARGET,SAVE        ::stdpara
type,extends(typgalaxy_para)::typspiral
       DOUBLE COMPLEX                   ::w
       INTEGER                          ::mode
       INTEGER                          ::N
       DOUBLE COMPLEX,ALLOCATABLE       ::u(:,:) 
       DOUBLE COMPLEX,ALLOCATABLE       ::h1(:)
       DOUBLE COMPLEX,ALLOCATABLE       ::phi1r(:)
       DOUBLE COMPLEX,ALLOCATABLE       ::k3(:)
       DOUBLE PRECISION,ALLOCATABLE     ::r(:)
       DOUBLE PRECISION                 ::rmax
       DOUBLE PRECISION                 ::rmin
       DOUBLE PRECISION                 ::fortoone
       DOUBLE COMPLEX                   ::error
       LOGICAL                          ::ucaled    = .false.
       LOGICAL                          ::h1caled   = .false.
       LOGICAL                          ::phi1rcaled= .false.
       LOGICAL                          ::bndu0     = .true.
       LOGICAL                          ::winit     = .false.
       CONTAINS
       PROCEDURE,NOPASS::init           =>spiral_init
       PROCEDURE::printu                =>spiral_printu
       PROCEDURE::printh1               =>spiral_printh1
       PROCEDURE::printr                =>spiral_printr
       PROCEDURE::final                 =>spiral_final
       PROCEDURE::readw                 =>spiral_readw
       PROCEDURE::printk3               =>spiral_printk3
endtype
CONTAINS

SUBROUTINE para_print(this)
IMPLICIT NONE
class(typgalaxy_para),intent(in)           ::this
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
namelist /paralist/ Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2
if(.not.this.inited)then
        write(0,*)'para not init,print failed'
else
      Lh        = this.para(1)
      rhoh      = this.para(2)
      Mb        = this.para(3)
      rb        = this.para(4)
      dM        = this.para(5)
      da        = this.para(6)
      db        = this.para(7)
      Qod       = this.para(8)
      q         = this.para(9)
      rq        = this.para(10)
      a1        = this.para(11)
      a2        = this.para(12)
      M1        = this.para(13)
      M2        = this.para(14)
      write(6,*)'para print:'
      write(6,nml=paralist)
endif
ENDSUBROUTINE

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
!NAME LIST
namelist /paralist/ Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2

!$OMP CRITICAL
open(10,file='para.list')
read(10,nml=paralist)
close(10)
!$OMP END CRITICAL

stdpara.para = (/Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2/)
stdpara.inited = .true.
ENDSUBROUTINE

SUBROUTINE cpfromstd(this)
IMPLICIT NONE
type(typgalaxy_para)                    ::this
this = stdpara
ENDSUBROUTINE

SUBROUTINE spiral_readw(this,mode)
IMPLICIT NONE
class(typspiral)                         ::this
!pspd from readin
DOUBLE PRECISION                        ::w(4)
INTEGER                                 ::mode
LOGICAL                                 ::bndu0(2)
namelist /spiralnml/ w,bndu0

this.mode = mode
!$OMP CRITICAL
open(10,file='para.list')
read(10,nml=spiralnml)
close(10)
!$OMP END CRITICAL
!this.bndu0 = bndu0(mode)

this.w = dcmplx(w(mode*2-1),w(mode*2))
this.winit = .true.
ENDSUBROUTINE

SUBROUTINE spiral_init(this,n,domain,para,mode)
USE OMP_LIB
IMPLICIT NONE
type(typspiral)                         ::this
type(typgalaxy_para)                    ::para
INTEGER                                 ::n,mode
DOUBLE PRECISION                        ::domain
if(.not.allocated(this.u))ALLOCATE(this.u(3,4*n))
if(.not.allocated(this.h1))ALLOCATE(this.h1(4*n))
if(.not.allocated(this.phi1r))ALLOCATE(this.phi1r(2*n))
if(.not.allocated(this.r))ALLOCATE(this.r(4*n))
if(.not.allocated(this.k3))ALLOCATE(this.k3(4*n))
!ALLOCATE(this.u(3,4*n))
!ALLOCATE(this.h1(4*n))
!ALLOCATE(this.phi1r(2*n))
!ALLOCATE(this.r(4*n))
this.rmin = 0.d0 
this.rmax = 1.5d0*domain
this.n    = 4*n
this.para = para.para
this.inited = .true.
ENDSUBROUTINE

SUBROUTINE spiral_final(this)
IMPLICIT NONE
class(typspiral)                         ::this
DEALLOCATE(this.u)
DEALLOCATE(this.h1)
DEALLOCATE(this.phi1r)
DEALLOCATE(this.r)
DEALLOCATE(this.k3)
ENDSUBROUTINE

SUBROUTINE spiral_printu(this)
IMPLICIT NONE
class(typspiral),INTENT(IN)             ::this
INTEGER                                 ::i
DO i = 1, this.n
        write(6,*)real(this.u(1,i)),abs(this.u(2,i))
ENDDO
ENDSUBROUTINE

SUBROUTINE spiral_printh1(this)
IMPLICIT NONE
class(typspiral),INTENT(IN)             ::this
INTEGER                                 ::i
DO i = 1, this.n
        write(6,*)this.r(i),abs(this.h1(i))
ENDDO
ENDSUBROUTINE

SUBROUTINE spiral_printr(this)
IMPLICIT NONE
class(typspiral),INTENT(IN)             ::this
INTEGER                                 ::i
DO i = 1, this.n
        write(6,*)this.r(i)
ENDDO
ENDSUBROUTINE

SUBROUTINE spiral_printk3(this)
IMPLICIT NONE
class(typspiral),INTENT(IN)             ::this
INTEGER                                 ::i
DO i = 1, this.n
        write(6,*)this.r(i),this.k3(i)
ENDDO
ENDSUBROUTINE

ENDMODULE

MODULE STELLARDISK
USE STELLARDISK_MODEL,only:stdpara,typspiral
USE PLOTTING
IMPLICIT NONE
DOUBLE PRECISION,PARAMETER::GravConst   = 4.3d-6 
DOUBLE PRECISION,PARAMETER::g           = 4.3d0
DOUBLE PRECISION,PARAMETER::pi          = 4.d0*atan(1.d0)
CONTAINS

SUBROUTINE FindSpiral(spiral)
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral

if(.not.spiral.inited)then
        write(0,*)'spiral object not initialized, stop'
        stop
endif
if(.not.spiral.winit)then
        write(0,*)'pspd not initialized,stop!!'
        stop
endif
!find EigenFunction
CALL findu(spiral)
!!find h1
CALL findh1(spiral)
!!Find phi1 along r
!CALL FindPhi1
!!Save calculation results
ENDSUBROUTINE
 
SUBROUTINE k3sqrtlog(spiral)
USE STELLARDISK_MODEL,only:stdpara,typspiral
IMPLICIT NONE                           
type(typspiral)                         ::spiral
INTEGER                                 ::i
DOUBLE PRECISION                        ::r
CHARACTER(len=32)                       ::ch
write(ch,"(A5,I1,A4)")'r-dep',spiral.mode,'.dat'
open(10,file=ch)
DO i = 2, spiral.n,2
        r = spiral.r(i)
        write(10,'(6(1XE15.6))')spiral.r(i),abs(spiral.u(2,i)),abs(spiral.h1(i))/snsd(r,spiral)**2*sigma0(r,spiral),abs(spiral.h1(i))
        !r, u, sigma1,h1
enddo
close(10)

!CALL rlog

ENDSUBROUTINE

!SUBROUTINE rlog
!USE plotting,only:plotlog
!IMPLICIT NONE
!DOUBLE PRECISION,ALLOCATABLE            ::dat(:,:)
!DOUBLE PRECISION                        ::r
!INTEGER                                 ::i,m
!
!!r, rotation, k3, Q
!m = 4
!
!ALLOCATE(dat(4,spiral.n))
!dat(1,:) = spiral.r
!do i = 1, spiral.n
!        r=spiral.r(i)
!        dat(2,i) = Omega(r)
!        dat(3,i) = dble(k3sqrt(r))
!        dat(4,i) = ToomreQ(r)
!enddo
!
!CALL plotlog(dat,m,spiral.n)
!DEALLOCATE(dat)
!ENDSUBROUTINE

function ToomreQ(r,spiral)
DOUBLE PRECISION                        ::Q,r,Qod,ToomreQ,rq
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)
para=>spiral.para
Qod = para(8)
q   = para(9)
rq  = para(10)
ToomreQ = Qod*(1.d0 + q*dexp(-r**2/rq**2))
endfunction

function nu(r,spiral)
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE COMPLEX                          ::nu   
DOUBLE PRECISION                        ::r
DOUBLE PRECISION                        ::m 
m = 2.d0
nu = (spiral.w-m*Omega(r,spiral))/kappa(r,spiral)
endfunction
 
function k3sqrt(r,spiral)
USE STELLARDISK_MODEL
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE COMPLEX            k3sqrt
DOUBLE PRECISION,INTENT(in)::r
DOUBLE PRECISION          ::rr


k3sqrt = (dcmplx(kappa(r,spiral)/snsd(r,spiral)))**2*(dcmplx(ToomreQ(r,spiral))**-2 &
       - 1.d0 + nu(r,spiral)**2 + 0.25d0*curF(r,spiral)**2*ToomreQ(r,spiral)**2)

!k3sqrt = (dcmplx(kappa(r)/snsd(r)))**2*(dcmplx(ToomreQ(r))**-2  &
!         - 1.d0 + nu(r)**2)

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

function snsd(r,spiral)
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION                        ::r,snsd
snsd = ToomreQ(r,spiral)*pi*GravConst*sigma0(r,spiral)/kappa(r,spiral)
ENDFUNCTION
 
function Sigma0(rr,spiral)
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)
DOUBLE PRECISION,INTENT(IN)::rr
DOUBLE PRECISION          ::Sigma0,r
DOUBLE PRECISION          ::m,a,b
DOUBLE PRECISION          ::BOUND,EPSREL,EPSABS
DOUBLE PRECISION          ::ans
DOUBLE PRECISION          ::ABSERR
INTEGER                   ::NEVAL,IERR,LIMIT,LENW,LAST,INF
DOUBLE PRECISION,ALLOCATABLE ::WORK(:)
INTEGER,ALLOCATABLE       ::IWORK(:)

para=>spiral.para

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

10 sigma0 = LauDiskSigma0(r,spiral)

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

FUNCTION LauDiskSigma0(r,spiral)
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION                ::LauDiskSigma0
DOUBLE PRECISION                ::r
!Lau Disk
DOUBLE PRECISION                ::x1,x2
DOUBLE PRECISION                ::a1,a2,M1,M2
DOUBLE PRECISION,POINTER                ::para(:)
para=>spiral.para
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

function kappa(r,spiral)
USE NUM
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION  kappa
DOUBLE PRECISION,INTENT(IN)             ::r
DOUBLE PRECISION  dr
DOUBLE PRECISION  dOmega
!if(r.lt.zerolimit)then
!        kappa = 2.d0*Omega(r)
!        return
!endif
kappa = sqrt(4.d0*Omega(r,spiral)**2*(1.d0+r/(2.d0*Omega(r,spiral))*dfunc(Omega,r,spiral)))
if(isnan(kappa))then
        print *,'kappa exception catched'
        print *,'r,Omega,dfunc(Omega,r):'
        print *,r,Omega(r,spiral),dfunc(Omega,r,spiral)
        CALL XERMSG('k3sqrt','kappa','kappa is nan.',-97,2)
endif       
ENDFUNCTION

FUNCTION Omega(r,spiral)
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
!spiral
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)=>null()

para=>spiral.para

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
        !!ignore by nightly test
!       Omega = 4.d0/3.d0*pi*GravConst*(MB+RHOH) &
!             + 16.d0/105*GravConst*(M1/a1**3-M2/a2**3)*120.d0
        Omega = Mb*GravConst/rb**3 &
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
goto 200
Mb   = para(3)
rb   = para(4)
gBulge = (-r/sqrt(1.d0+r**2/rb**2)/rb+dasinh(r/rb))/r**2

gBulge = 4.d0*pi*rb**3*Mb*gBulge*GravConst
VBulge = sqrt(r*gBulge)

!Lau's Bulge for nightly test
200 VBulge = sqrt(MB*GravConst/rb**3*(1.d0+r**2/rb**2)**-1.5)*r


!!Disk
!dM     = para(5)
!da     = para(6)
!db     = para(7)
!!VDisk  = sqrt(dfunc(pDisk,r)*r)
!VDisk  = sqrt(dpDisk(r)*r)

!LauDisk 
VDisk  = VLauDisk(r)
!!ignore by nightly test
!Omega = sqrt(VHalo**2+VBulge**2+dble(VDisk**2))/r
Omega = sqrt(VBulge**2+dble(VDisk**2))/r
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

function dfunc(func,r,spiral)
USE STELLARDISK_MODEL
!
! Forward differential
!
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE PRECISION                ::dr
DOUBLE PRECISION,PARAMETER      ::coe(3)=(/-1.5d0,2.d0,-0.5d0/)
DOUBLE PRECISION                ::dfunc,ans,funcs(3)
INTEGER                         ::i
interface 
        function func(x,spiral)
        USE STELLARDISK_MODEL
        DOUBLE PRECISION        ::func
        DOUBLE PRECISION        ::x
        type(typspiral),TARGET                  ::spiral
        ENDFUNCTION func
endinterface        

dr = epsilon(r)**0.3*max(r,epsilon(0d0))
funcs(1) = func(r,spiral)
funcs(2) = func(r+dr,spiral)
funcs(3) = func(r+2.d0*dr,spiral)
ans = dot_product(funcs,coe)/dr
dfunc = ans

!ans = 0.d0
!ans =        -3.d0/2.d0*func(r)
!ans = ans +        2.d0*func(r+dr)
!ans = ans +  -1.d0/2.d0*func(r+2.d0*dr)
!dfunc = ans/dr
if(isnan(dfunc))CALL XERMSG('k3sqrt','dfunc','dfunc is nan.',-94,0)
endfunction
 
function curf(r,spiral)
USE NUM
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION                ::curf
DOUBLE PRECISION                ::r,tmp,s
INTEGER                         ::m=2

if(r.lt.zerolimit)then
        curf = &
        2.d0*dble(m)*(pi*GravConst*sigma0(0.d0,spiral))/kappa(0.d0,spiral)**2*sqrt(-2.d0*Omega2(spiral)/(Omega(0.d0,spiral)))
        return
elseif(r.lt.3.d-2)then
        curf = &
        2.d0*dble(m)*(pi*GravConst*sigma0(r,spiral))/kappa(r,spiral)**2 &
        *sqrt(-(2.d0*r*Omega2(spiral))/(r*Omega(r,spiral)))
else
        s = -r/Omega(r,spiral)*dfunc(Omega,r,spiral)
!       curf = &
!       2.d0*dble(m)*(pi*GravConst*sigma0(r,spiral))/kappa(r,spiral)**2*sqrt(-2.d0*dfunc(Omega,r,spiral)/(2.d0*r*Omega(r,spiral)+dfunc(Omega,r,spiral)))
        curf = &
        2.d0*dble(m)*(pi*GravConst*sigma0(r,spiral))/kappa(r,spiral)**2/r/sqrt(1.d0/s-0.5d0)
endif

CALL CheckResult
CONTAINS
FUNCTION Omega2(spiral)
USE STELLARDISK_MODEL
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION          ::Omega2
!Halo
DOUBLE PRECISION          ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION          ::rb,Mb,gBulge,VBulge
!disk
DOUBLE PRECISION          ::dM,da,db
DOUBLE COMPLEX            ::VDisk
DOUBLE PRECISION          ::a1,a2,M1,M2
!para
DOUBLE PRECISION,POINTER                ::para(:)
para=>spiral.para

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

!ignored by nightly test
!Omega2 = 4.d0*pi*GravConst*rhoh/5.d0/Lh**2 + 6.d0/5.d0*pi*GravConst*Mb/rb**2 &
!       + 8.d0*GravConst/105.d0*(M1/a1**5-M2/a2**5)*1080.d0
Omega2 = 3.d0/2.d0*GravConst*Mb/rb**5 &
       + 8.d0*GravConst/105.d0*(M1/a1**5-M2/a2**5)*1080.d0
Omega2 = -Omega2/2.d0/Omega(0.d0,spiral)
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
        print *,r,Omega(r,spiral),dfunc(Omega,r,spiral)
        CALL XERMSG('k3sqrt','curf',errormsg,-96,2)
endif
ENDSUBROUTINE
ENDFUNCTION
 
SUBROUTINE findu(spiral)
USE RK,only:rk4
USE STELLARDISK_MODEL
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE COMPLEX                  ::ui(3)
DOUBLE PRECISION                ::a,b
INTEGER                         ::n,i

spiral.u = (1.d0,0.d0)*0.d0
a = spiral.rmin
b = spiral.rmax
if(spiral.bndu0)then
        ui = (/dcmplx(a),dcmplx(0.d0),2.d0*sqrt(-q(0.d0))/)
else
        ui = (/dcmplx(a),dcmplx(1.d0,0.d0),dcmplx(0.d0)/)
endif
CALL rk4(a,b,spiral.n,p,q,p,spiral.u,ui)
spiral.ucaled = .true.
spiral.r(:) = real(spiral.u(1,:))
spiral.error = error(spiral)
DO i = 1, spiral.n
        spiral.k3(i) = k3sqrt(spiral.r(i),spiral)
ENDDO
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
q = k3sqrt(r,spiral)
!q =  (1)
!q = dcmplx(1.d0,r**2)
ENDFUNCTION

ENDSUBROUTINE

SUBROUTINE findh1(spiral)
IMPLICIT NONE
type(typspiral)                         ::spiral
DOUBLE COMPLEX                  ::u
DOUBLE COMPLEX                  ::rad
DOUBLE PRECISION                ::h,r,expp
INTEGER                         ::i
do i = 1, spiral.n
        r   = spiral.r(i)
        u   = spiral.u(2,i)
        rad = sqrt(kappa(r,spiral)**2*(1.d0-nu(r,spiral)**2)/sigma0(r,spiral)/r)
        spiral.h1(i) = u*rad*exp(-0.5d0*(0.d0,1.d0)*ExpPart(r))
enddo

!CALL refineh1
spiral.h1caled = .true.

CONTAINS

Function ExpPart(r) result(ans)
USE NUM
        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(IN)     ::r
        DOUBLE PRECISION                ::a,rr
        DOUBLE PRECISION                ::ans
        DOUBLE PRECISION                ::ferr = 1.d-15
        INTEGER                         ::IERR,K=6000
        a = zerolimit
        rr = r
        IERR = 0
        if(rr.eq.0.d0)then
                ans = 0.d0
        else
                CALL DGAUS8(F,a,rr,fERR,ans,IERR)
!               CALL DQNC79(F,a,rr,fERR,ans,IERR,K)
        endif

ENDFUNCTION ExpPart

FUNCTION F(r)
DOUBLE PRECISION                ::F,r
F = Sigma(r,spiral)
ENDFUNCTION


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

SUBROUTINE FindPhi1(spiral)
USE STELLARDISK_MODEL
IMPLICIT NONE
DOUBLE COMPLEX                  ::k(4)
DOUBLE PRECISION                ::r,h
INTEGER                         ::i,j,l,n
type(typspiral),TARGET                  ::spiral

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
!
!dsimplifiedPoisson = -phi/(2.d0*r)+(0.d0,1.d0)*cmplx(Sigma(r),0)*h
!dsimplifiedPoisson = dsimplifiedPoisson + 3.75d0/(0.d0,1.d0)/Sigma(r)/r**2*phi
!
!!!test case 
!!!dsimplifiedPoisson =  -phi*r + (0.d0,1.d0)*r
!!!bnd condition is phi = 1+i at r=0
!!!solution is
!!!phi = exp(-r**2/2)+ i
!!!dsimplifiedPoisson = -phi*r
ENDFUNCTION

function Sigma(r,spiral) RESULT(ans)
USE STELLARDISK_MODEL
!This is NOT related to density
IMPLICIT NONE
type(typspiral)                         ::spiral
DOUBLE PRECISION                ::ans
DOUBLE PRECISION,INTENT(IN)     ::r
ans = 2.d0*pi*GravConst*sigma0(r,spiral)/snsd(r,spiral)**2
ENDFUNCTION

SUBROUTINE FindForce(force,r,th)
USE STELLARDISK_MODEL
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
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

FUNCTION sigma1(r,th,spiral)
USE STELLARDISK_MODEL
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
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
sigma1 = real(hh1*sigma0(r,spiral)/snsd(r,spiral)**2*exp(-2.d0*th*(0.d0,1.d0)))


ENDFUNCTION

FUNCTION phi1(r,th)
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
DOUBLE PRECISION                ::phi1
DOUBLE PRECISION,INTENT(IN)     ::r,th
DOUBLE COMPLEX                  ::uu,hh1
DOUBLE PRECISION                ::rad
INTEGER                         ::i,j,k,l,n
!!interploting u at non-grid point r
!n = spiral.n
!do i = 1, n*2
!        if(spiral.r(2*i).gt.r)then
!                exit
!        endif
!enddo
!!uu = (r-u(1,i-1))/(u(1,i)-u(1,i-1))*(u(2,i)-u(2,i-1)) + u(2,i-1)
!hh1 =(r-spiral.u(1,2*i-2))/(spiral.u(1,2*i)-spiral.u(1,2*i-2))*(spiral.phi1r(i)-spiral.phi1r(i-1)) + spiral.phi1r(i-1)
!
!
!!find potential
!phi1 = real(hh1*exp(-2.d0*th*(0.d0,1.d0)))
!
!
ENDFUNCTION

SUBROUTINE ENDSTELLARDISK
!DEALLOCATE(spiral.u)
!DEALLOCATE(spiral.h1)
!DEALLOCATE(spiral.phi1r)
!DEALLOCATE(spiral.r)
ENDSUBROUTINE

function error(spiral)
USE STELLARDISK_MODEL
IMPLICIT NONE          
type(typspiral)                         ::spiral
DOUBLE COMPLEX          ::error
DOUBLE COMPLEX          ::uu(3)
DOUBLE PRECISION        ::h=10d-5   ,r
DOUBLE PRECISION        ::RE,AE,B,C,RR
INTEGER                 ::l,IFLAG

B = 6.d0
C = 20.d0
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
error  = -(0.d0,1.d0)*sqrt(k3sqrt(r,spiral))
error  = error -                           &
0.5d0/sqrt(k3sqrt(r,spiral))*               &                
(sqrt(k3sqrt(r+h,spiral))-sqrt(k3sqrt(r-h,spiral)))/(2.d0*h)
error = uu(3)/uu(2) -error
CONTAINS
function four21(r)
IMPLICIT NONE
DOUBLE PRECISION                ::four21,r
four21 = (real(spiral.w)- 2.d0*Omega(r,spiral))/kappa(r,spiral) - 0.5d0
ENDFUNCTION
        
endfunction

ENDMODULE STELLARDISK
