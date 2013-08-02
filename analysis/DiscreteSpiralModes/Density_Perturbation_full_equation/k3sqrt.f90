MODULE NUM
DOUBLE PRECISION,PARAMETER                      ::zerolimit=1.d-6
ENDMODULE NUM

MODULE STELLARDISK_MODEL
!!Define Types
type   typgalaxy_para
       DOUBLE PRECISION                 ::para(14)
       LOGICAL                          ::inited=.false.
       CONTAINS 
       PROCEDURE                        ::print         =>para_print
       PROCEDURE                        ::readstd       =>READINSTDPARA
endtype
type,extends(typgalaxy_para)            ::typspiral
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
       DOUBLE PRECISION                 ::co              ! position of corotation
       DOUBLE PRECISION                 ::fortoone        ! position of outer four to one
       DOUBLE COMPLEX                   ::error
       DOUBLE PRECISION                 ::phase     = 0.d0! starting angle for 2D plot
       DOUBLE PRECISION                 ::dr              ! deparation btwn two points
       DOUBLE PRECISION                 ::pspd            ! pattern speed
       LOGICAL                          ::ucaled    = .false.
       LOGICAL                          ::cocaled   = .false.
       LOGICAL                          ::h1caled   = .false.
       LOGICAL                          ::bndu0     = .true.
       LOGICAL                          ::winit     = .false.
       CONTAINS
       PROCEDURE,PASS                   ::init          =>spiral_init
       PROCEDURE,PASS                   ::free          =>spiral_final
       PROCEDURE,PASS                   ::readw         =>spiral_readw
       PROCEDURE,PASS                   ::setw          =>spiral_setw
endtype
!!STELLAR VARIABLES
type(typgalaxy_para),TARGET,SAVE        ::stdpara
CHARACTER(*),PARAMETER                  ::parafnm = 'para.list'
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

SUBROUTINE READINSTDPARA(this)
IMPLICIT NONE
CLASS(typgalaxy_para)           ::this
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
open(10,file=parafnm)
read(10,nml=spiralnml)
close(10)
!$OMP END CRITICAL
!this.bndu0 = bndu0(mode)

this.w = dcmplx(w(mode*2-1),w(mode*2))
this.winit = .true.
this.pspd = real(this.w)/2.d0
ENDSUBROUTINE
 
SUBROUTINE spiral_setw(this,w,mode)
IMPLICIT NONE
class(typspiral)                        ::this
DOUBLE COMPLEX                          ::w
INTEGER                                 ::mode

this.mode = mode
this.w = w 
this.winit = .true.
this.pspd = real(this.w)/2.d0
ENDSUBROUTINE

SUBROUTINE spiral_init(this,n,domain,para,mode)
IMPLICIT NONE
class(typspiral)                        ::this
type(typgalaxy_para)                    ::para
INTEGER                                 ::n,mode
DOUBLE PRECISION                        ::domain
IF(.NOT.ALLOCATED(this.u))ALLOCATE(this.u(3,2*n))
IF(.NOT.ALLOCATED(this.h1))ALLOCATE(this.h1(2*n))
IF(.NOT.ALLOCATED(this.r))ALLOCATE(this.r(2*n))
IF(.NOT.ALLOCATED(this.k3))ALLOCATE(this.k3(2*n))
this.rmin = 0.d0 
this.rmax = 1.5d0*domain
this.n    = 2*n
this.para = para.para
this.inited = .true.
this.dr   = (this.rmax - this.rmin)/dble(this.n)
this.co   = -100.d0
ENDSUBROUTINE

SUBROUTINE spiral_final(this)
IMPLICIT NONE
class(typspiral)                         ::this
IF(ALLOCATED(this.u))DEALLOCATE(this.u)
IF(ALLOCATED(this.h1))DEALLOCATE(this.h1)
IF(ALLOCATED(this.r))DEALLOCATE(this.r)
IF(ALLOCATED(this.k3))DEALLOCATE(this.k3)
ENDSUBROUTINE

ENDMODULE

MODULE STELLARDISK
USE STELLARDISK_MODEL,only:stdpara,typspiral
IMPLICIT NONE
DOUBLE PRECISION,PARAMETER::GravConst   = 4.3d-6 
DOUBLE PRECISION,PARAMETER::g           = 4.3d0
DOUBLE PRECISION,PARAMETER::pi          = 4.d0*atan(1.d0)
CONTAINS

SUBROUTINE FindSpiral(spiral)
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral

if(.not.spiral.inited)then
        write(0,*)'spiral object not initialized, stop. init flag:',spiral.inited
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
ENDSUBROUTINE
 
FUNCTION ToomreQ(r,spiral)
DOUBLE PRECISION                        ::Q,r,Qod,ToomreQ,rq
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)
DOUBLE PRECISION                        ::f
para=>spiral.para
Qod = para(8)
q   = para(9)
rq  = para(10)

ToomreQ = Qod*(1.d0 + q*dexp(-r**2/rq**2) + 1.2d0*dexp(-r**2/0.8**2)) 

!IF(r<rq)THEN
!        ToomreQ = Cos(r*pi/rq)*q+Qod
!        ToomreQ = ToomreQ + 1.0d0*exp(-r**2/0.8d0)
!ELSE
!        ToomreQ = Qod - q
!ENDIF
ENDFUNCTION

FUNCTION nu(r,spiral)
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE COMPLEX                          ::nu   
DOUBLE PRECISION                        ::r
DOUBLE PRECISION,PARAMETER              ::m =2.d0
nu = (spiral.w-m*Omega(r,spiral))/kappa(r,spiral)
ENDFUNCTION
 
function k3sqrt(r,spiral,outk)
USE STELLARDISK_MODEL
USE NUM
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE COMPLEX            ::k3sqrt,k(6)
DOUBLE COMPLEX,OPTIONAL   ::outk(6)
DOUBLE COMPLEX            ::w
DOUBLE PRECISION,INTENT(in)::r
DOUBLE PRECISION          ::rr
DOUBLE PRECISION,PARAMETER::m = 2.d0
DOUBLE PRECISION          ::Omega_, ToomreQ_, snsd_, sigma_, kappa_
DOUBLE PRECISION          ::sigma0_,sigma02_,sigma00_,sigma04_
DOUBLE PRECISION          ::kappa2_,kappa0_,Omega0_,Omega2_,Omega4_
DOUBLE COMPLEX            ::nu_
LOGICAL                   ::MASK(6)=.true.!this is for sum of part of k
INTEGER                   ::I

!!set co
CALL set_co(spiral)

Omega_   = Omega(r,spiral)
ToomreQ_ = ToomreQ(r,spiral)
snsd_    = snsd(r,spiral)
Sigma0_  = Sigma0(r,spiral)
Sigma_   = Sigma(r,spiral)
kappa_   = kappa(r,spiral)
nu_      = nu(r,spiral)
w = spiral.w
Omega0_  = OmegaCoefficent(spiral,0)
Omega2_  = OmegaCoefficent(spiral,2)
Omega4_  = OmegaCoefficent(spiral,4)
kappa0_  = kappa0(spiral)
kappa2_  = kappa2(spiral)
sigma00_ = Sigma0(0.d0,spiral)
sigma02_ = LauDiskSigma02(spiral)
sigma04_ = LauDiskSigma04(spiral)

k(:)     = 0.d0
k3sqrt   = 0.d0

k(1)     = (dcmplx(kappa_/snsd_))**2*(dcmplx(ToomreQ_)**-2 &
           - 1.d0 + nu_**2 + 0.25d0*curF(r,spiral)**2*ToomreQ_**2)

IF(r.le.2.0d-1)THEN
        !!r too small
        !!k(2) hydro
        !k(2) = (-2.d0*sigma02_/sigma00_+4d0*Omega2_/Omega0_) &
        !     + (2d0*sigma02_**2/sigma00_**2 - 4d0*sigma04_/sigma00_+4d0*(-2d0*Omega2_**2+3d0*Omega0_*Omega4_)/Omega0_**2)*r**2
        !k(2) = -2.d0*m*Omega_/(spiral.w-m*Omega_)*k(2)
        k(3) = - dcmplx(0.d0,1.d0)*Sigma_*5d-1/f2(r,spiral)*dcfunc(f2,r,spiral)

        k(4) = sigma02_**2*w**2*(w-4d0*Omega0_)**2 - 8d0*sigma00_*sigma02_*w*(w-4d0*Omega0_)*(w+Omega0_)*Omega2_ &
             + 4d0*sigma00_*(-sigma04_*w**2*(w-4d0*Omega0_)**2-4d0*sigma00_*(4d0*w**2+2d0*w*Omega0_+3d0*Omega0_**2)*Omega2_**2- &
             4d0*sigma00_*w*(w-4d0*Omega0_)*(w+2d0*Omega0_)*Omega4_)
        k(4) = k(4)/sigma00_**2/w**2/(w-4d0*Omega0_)**2*r**2

        k(5) = (3d0*Omega2_ + m*nu_*Omega2_)*5d-1/Omega0_ &
             + (-13d0*Omega2_**2-4d0*m*nu_*Omega2_**2+16d0*Omega0_*Omega4_+4d0*m*nu_*Omega0_*Omega4_)*0.25d0/Omega0_**2*r**2
        k(5) = -k(5)*(4d0*m*nu_/(1.d0-nu_**2))
        !! plasma 
        k(6) = (-2.d0*sigma02_/sigma00_+4d0*Omega2_/Omega0_) &
             + (2d0*sigma02_**2/sigma00_**2 - 4d0*sigma04_/sigma00_+4d0*(-2d0*Omega2_**2+3d0*Omega0_*Omega4_)/Omega0_**2)*r**2
        k(6) =  2.d0*m*Omega_/kappa_*b(r,spiral)*zbar(b(r,spiral)*nu_)*k(6)
ELSE
        !!hydro 
        !k(2) = - 2.d0*Omega_*m/(spiral.w-m*Omega_)/r*dcmplx(dfunc(f1,r,spiral)/f1(r,spiral))
        k(3) = - dcmplx(0.d0,1.d0)*Sigma_*5d-1/f2(r,spiral)*dcfunc(f2,r,spiral)
        k(4) = + 3.d0/4.d0/r**2  - cd2func(f3,r,spiral)/f3(r,spiral)  
        k(5) = - 4.d0*m*Omega_/kappa_**2*nu_*(m*nu_*dfunc(Omega,r,spiral)+dfunc(kappa,r,spiral))/r/(1.d0-nu_**2)
        !! plasma 
        k(6) =  2.d0*Omega_*m/r/kappa_*b(r,spiral)*zbar(b(r,spiral)*nu_)*dcmplx(dfunc(f1,r,spiral)/f1(r,spiral))
ENDIF

IF(PRESENT(outk))outk = k

MASK(2)= .false.
k3sqrt = sum(k,MASK)
!write(*,'(7(D15.5))')r,real(k3sqrt),real(k(1)),(real(k(i)),I = 3,6)
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

FUNCTION f1(r,spiral)
type(typspiral)                         ::spiral
DOUBLE PRECISION        ::f1
DOUBLE PRECISION        ::r

f1 = kappa(r,spiral)**2/sigma0(r,spiral)/Omega(r,spiral)

ENDFUNCTION

FUNCTION f2(r,spiral)
type(typspiral)                         ::spiral
DOUBLE PRECISION        ::f2
DOUBLE PRECISION        ::r

f2 = ToomreQ(r,spiral)**2*(1.d0-nu(r,spiral)**2)

ENDFUNCTION

FUNCTION f3(r,spiral)
type(typspiral)                         ::spiral
DOUBLE COMPLEX          ::f3
DOUBLE PRECISION        ::r,xr,xi

f3 = sqrt(sigma0(r,spiral)/r/kappa(r,spiral)**2/(dcmplx(1.d0)-nu(r,spiral)**2))

ENDFUNCTION

FUNCTION f4(r,spiral)
type(typspiral)                         ::spiral
DOUBLE PRECISION        ::f4
DOUBLE PRECISION        ::r

f4 = Omega(r,spiral)

ENDFUNCTION

FUNCTION zbar(z_in)
IMPLICIT NONE
DOUBLE COMPLEX                          ::zbar,z_in,z
DOUBLE PRECISION                        ::a,b
LOGICAL                                 ::FLAG

z = -z_in
CALL wofz(real(z),imag(z),a,b,FLAG)
zbar = -dcmplx(a,b)*dcmplx(0d0,sqrt(pi))

ENDFUNCTION

FUNCTION b(r,spiral)
IMPLICIT NONE
type(typspiral)                         ::spiral
DOUBLE PRECISION                        ::r
DOUBLE COMPLEX                          ::b,invdnu

IF(r.lt.0.2d0)THEN
        invdnu = -2d0*Omega0_**2/(3d0*Omega2_*w-2d0*Omega2_*Omega0_) &
             + Omega0_*((-19d0*w+14d0*Omega0_)*Omega2_**2+16d0*(w-Omega0_)*Omega0_*Omega4_)/Omega2_**2/(3d0*w-2d0*Omega0_)**2*r**2
             b = kappa_/sqrt(2d0)/spiral.co/snsd_*invdnu
ELSE
        b = kappa_*r/sqrt(2d0)/spiral.co/snsd_/cdfunc(nu,r,spiral)
ENDIF
ENDFUNCTION

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

Omega2 = 4.d0*pi*GravConst*rhoh/5.d0/Lh**2 + 6.d0/5.d0*pi*GravConst*Mb/rb**2 &
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

function snsd(r,spiral)
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION                        ::r,snsd
snsd = ToomreQ(r,spiral)*pi*GravConst*sigma0(r,spiral)/kappa(r,spiral)
ENDFUNCTION
 
function Sigma0(rr,spiral)
!This is surface density of disk
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

FUNCTION LauDiskSigma02(spiral)
!!coeficient of r**2 in disk density
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION                        ::LauDiskSigma02
!Lau Disk
DOUBLE PRECISION                        ::a1,a2,M1,M2
DOUBLE PRECISION,POINTER                ::para(:)
para=>spiral.para
a1 = para(11)
a2 = para(12)
M1 = para(13)
M2 = para(14)

LauDiskSigma02 = (-M1/a1**4/2.d0+M2/a2**4/2.d0)*4.5d0/pi*11.d0

ENDFUNCTION

FUNCTION LauDiskSigma04(spiral)
!!coeficient of r**2 in disk density
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION                        ::LauDiskSigma04
!Lau Disk
DOUBLE PRECISION                        ::a1,a2,M1,M2
DOUBLE PRECISION,POINTER                ::para(:)
para=>spiral.para
a1 = para(11)
a2 = para(12)
M1 = para(13)
M2 = para(14)

LauDiskSigma04 = (M1/a1**6-M2/a2**6)*143d0/8d0*4.5d0/pi
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

function kappa0(spiral)
USE NUM
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION                        ::kappa0
kappa0 = 2*Omega(0d0,spiral)
ENDFUNCTION

function kappa2(spiral)
USE NUM
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION                        ::kappa2
kappa2 = 3d0*Omega2(spiral)
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

VLauDisk = (M1*GravConst/a1**3*H(x1)) - (M2*GravConst/a2**3*H(x2))
VLauDisk = sqrt(16.d0/105.d0*VLauDisk)
VLauDisk = VLauDisk*r
ENDFUNCTION

FUNCTION H(x)
IMPLICIT NONE
DOUBLE PRECISION                ::H,x
DOUBLE PRECISION,PARAMETER      ::C(5) = (/59.0625d0,26.25d0,16.875d0,11.25d0,6.5625d0/)
DOUBLE PRECISION                ::xs(5),xx
INTEGER                         ::i

xx    = x**-2
xs(1) = x**-3

DO i = 2, 5
        xs(i) = xs(i-1)*xx
ENDDO

H = dot_product(C,xs)

!H = 0.d0
!H = 59.0625d0/x**11 + H
!H = 26.25d0  /x**9  + H
!H = 16.875d0 /x**7  + H
!H = 11.25d0  /x**5  + H
!H = 6.5625d0 /x**3  + H

ENDFUNCTION

ENDFUNCTION

FUNCTION OmegaCoefficent(spiral,n)
IMPLICIT NONE
DOUBLE PRECISION                        ::OmegaCoefficent
type(typspiral),TARGET,INTENT(IN)       ::spiral
INTEGER,INTENT(IN)                      ::N

SELECT CASE(n)
        CASE(0)
                OmegaCoefficent = Omega0(spiral)
        CASE(2)
                OmegaCoefficent = Omega2(spiral)
        CASE(4)
                OmegaCoefficent = Omega4(spiral)
ENDSELECT
ENDFUNCTION

FUNCTION Omega0(spiral)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION          ::Omega0
type(typspiral),TARGET                  ::spiral
Omega0 = Omega(0d0,spiral)
ENDFUNCTION

FUNCTION Omega2(spiral)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION          ::Omega2
!Halo
DOUBLE PRECISION          ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION          ::rb,Mb,gBulge,VBulge
!disk
DOUBLE COMPLEX            ::VDisk
DOUBLE PRECISION          ::a1,a2,M1,M2
!spiral
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)=>null()
DOUBLE PRECISION                        ::tmp(3)

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

tmp = 0.d0
tmp(1) = 4.d0*pi*GravConst*rhoh/5.d0/Lh**2
tmp(2) = 6.d0*pi*GravConst*Mb/5.d0/rb**2
tmp(3) = 8.d0*GravConst/105d0*(M1/A1**5-M2/A2**5)*1080d0

Omega2 = sum(tmp)/Omega(0.d0,Spiral)/-2.d0
CONTAINS

SUBROUTINE CheckResult
IMPLICIT NONE
CHARACTER(len=72)                       ::errormsg
write(errormsg,"(E10.3)")0.d0
errormsg = trim(errormsg)
errormsg = 'Omega nan.@r = '//errormsg
if(isnan(Omega2))CALL XERMSG('k3sqrt','Omega2',errormsg,-99,2)
ENDSUBROUTINE


ENDFUNCTION

FUNCTION Omega4(spiral)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION          ::Omega4,Omega0_
!Halo
DOUBLE PRECISION          ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION          ::rb,Mb,gBulge,VBulge
!disk
DOUBLE COMPLEX            ::VDisk
DOUBLE PRECISION          ::a1,a2,M1,M2
!spiral
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)=>null()
DOUBLE PRECISION                        ::tmp(4)

para=>spiral.para

Omega0_ = Omega0(spiral)

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

tmp = 0.d0
tmp(1) = -Omega2(spiral)**2
tmp(2) = (4.d0*pi*GravConst*rhoh/7.d0/Lh**4)
tmp(3) = (15.d0*pi*GravConst*Mb/14.d0/rb**4)
tmp(4) = (2.d0*GravConst/105d0*(M1/A1**7-M2/A2**7)*(8280d0))

Omega4 = sum(tmp)*5d-1/Omega0_
CONTAINS

SUBROUTINE CheckResult
IMPLICIT NONE
CHARACTER(len=72)                       ::errormsg
write(errormsg,"(E10.3)")0.d0
errormsg = trim(errormsg)
errormsg = 'Omega nan.@r = '//errormsg
if(isnan(Omega4))CALL XERMSG('k3sqrt','Omega4',errormsg,-99,2)
ENDSUBROUTINE


ENDFUNCTION

FUNCTION Oomega(r,spiral)
!use summantion instead of put in the same subroutine
USE NUM
IMPLICIT NONE
DOUBLE PRECISION          ::Oomega
DOUBLE PRECISION,INTENT(IN)::r
!spiral
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)=>null()

para=>spiral.para

Oomega = bOmega(r,spiral)**2 &
       + dOmega(r,spiral)**2 &
       + hOmega(r,spiral)**2

Oomega = sqrt(Oomega)

ENDFUNCTION

FUNCTION bOmega(r,spiral)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION          ::bOmega
DOUBLE PRECISION,INTENT(IN)::r
!bulge
DOUBLE PRECISION          ::rb,Mb,gBulge,VBulge
!spiral
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)=>null()

para=>spiral.para

!Bulge
Mb   = para(3)
rb   = para(4)

!!Limit case when r->0:
if(r.lt.zerolimit)then
        bOmega = 4.d0/3.d0*pi*GravConst*(MB) 
        bOmega = bOmega**0.5
        return
endif


!Bulge
Mb   = para(3)
rb   = para(4)
gBulge = (-r/sqrt(1.d0+r**2/rb**2)/rb+dasinh(r/rb))/r**2

gBulge = 4.d0*pi*rb**3*Mb*gBulge*GravConst
VBulge = sqrt(r*gBulge)

bOmega = VBulge/r
!print *,r,V
!CALL CheckResult
CONTAINS 
SUBROUTINE CheckResult
IMPLICIT NONE
CHARACTER(len=72)                       ::errormsg
write(errormsg,"(E10.3)")r
errormsg = trim(errormsg)
errormsg = 'bOmega nan.@r = '//errormsg
if(isnan(bOmega))CALL XERMSG('k3sqrt','bOmega',errormsg,-99,2)
ENDSUBROUTINE

ENDFUNCTION

FUNCTION dOmega(r,spiral)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION          ::dOmega
DOUBLE PRECISION,INTENT(IN)::r
!disk
DOUBLE PRECISION          ::dM,da,db
DOUBLE COMPLEX            ::VDisk
DOUBLE PRECISION          ::a1,a2,M1,M2
!spiral
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)=>null()

para=>spiral.para

!Disk
a1 = para(11)
a2 = para(12)
M1 = para(13)
M2 = para(14)

!!Limit case when r->0:
if(r.lt.zerolimit)then
        dOmega = 16.d0/105*GravConst*(M1/a1**3-M2/a2**3)*120.d0
        dOmega = dOmega**0.5
        return
endif



!LauDisk 
VDisk  = VLauDisk(r)
dOmega = VDisk/r
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
if(isnan(dOmega))CALL XERMSG('k3sqrt','Omega',errormsg,-99,2)
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

FUNCTION hOmega(r,spiral)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION          ::hOmega
DOUBLE PRECISION,INTENT(IN)::r
!Halo
DOUBLE PRECISION          ::Lh,rhoh,gHalo,VHalo
!spiral
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,POINTER                ::para(:)=>null()

para=>spiral.para

!Halo
Lh   = para(1)
rhoh = para(2)
!!Limit case when r->0:
if(r.lt.zerolimit)then
        hOmega = 4.d0/3.d0*pi*GravConst*(RHOH) 
        hOmega = hOmega**0.5
        return
endif


!Halo
Lh   = para(1)
rhoh = para(2)
gHalo = 4.d0*Lh**2.d0*pi*rhoh*(r - Lh*atan(r/Lh))/r**2
gHalo = GravConst*gHalo
VHalo = sqrt(r*gHalo)


hOmega = VHalo/r
!print *,r,V
!CALL CheckResult
CONTAINS


SUBROUTINE CheckResult
IMPLICIT NONE
CHARACTER(len=72)                       ::errormsg
write(errormsg,"(E10.3)")r
errormsg = trim(errormsg)
errormsg = 'Omega nan.@r = '//errormsg
if(isnan(hOmega))CALL XERMSG('k3sqrt','Omega',errormsg,-99,2)
ENDSUBROUTINE


ENDFUNCTION

FUNCTION StellarOmega(r,spiral)
USE NUM
IMPLICIT NONE
DOUBLE PRECISION          ::StellarOmega
DOUBLE PRECISION,INTENT(IN)::r
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
        StellarOmega = 4.d0/3.d0*pi*GravConst*(MB) &
              + 16.d0/105*GravConst*(M1/a1**3-M2/a2**3)*120.d0
        StellarOmega = StellarOmega**0.5
        return
endif


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
StellarOmega = sqrt(VBulge**2+dble(VDisk**2))/r
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
if(isnan(StellarOmega))CALL XERMSG('k3sqrt','Omega',errormsg,-99,2)
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

function cdfunc(func,r,spiral)
USE STELLARDISK_MODEL
!
! Forward differential
!
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE PRECISION                ::dr
DOUBLE PRECISION,PARAMETER      ::coe(3)=(/-1.5d0,2.d0,-0.5d0/)
DOUBLE COMPLEX                  ::cdfunc,ans,funcs(3)
INTEGER                         ::i
interface 
        function func(x,spiral)
        USE STELLARDISK_MODEL
        DOUBLE COMPLEX          ::func
        DOUBLE PRECISION        ::x
        type(typspiral),TARGET                  ::spiral
        ENDFUNCTION func
endinterface        

dr = epsilon(r)**0.3*max(r,epsilon(0d0))
funcs(1) = func(r,spiral)
funcs(2) = func(r+dr,spiral)
funcs(3) = func(r+2.d0*dr,spiral)
ans = dot_product(coe,funcs)/dr
cdfunc = ans

!ans = 0.d0
!ans =        -3.d0/2.d0*func(r)
!ans = ans +        2.d0*func(r+dr)
!ans = ans +  -1.d0/2.d0*func(r+2.d0*dr)
!dfunc = ans/dr
if(isnan(real(cdfunc)))CALL XERMSG('k3sqrt','cdfunc','cdfunc is nan.',-94,0)
endfunction

function dcfunc(func,r,spiral)
USE STELLARDISK_MODEL
!
! Forward differential
!
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE PRECISION                ::dr
DOUBLE PRECISION,PARAMETER      ::coe(3)=(/-1.5d0,2.d0,-0.5d0/)
DOUBLE COMPLEX                  ::dcfunc,ans,funcs(3)
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
dr = 1d-3
funcs(1) = func(r,spiral)
funcs(2) = func(r+dr,spiral)
funcs(3) = func(r+2.d0*dr,spiral)
ans = dot_product(funcs,coe)/dr
dcfunc = ans

!ans = 0.d0
!ans =        -3.d0/2.d0*func(r)
!ans = ans +        2.d0*func(r+dr)
!ans = ans +  -1.d0/2.d0*func(r+2.d0*dr)
!dfunc = ans/dr
if(isnan(real(dcfunc)))CALL XERMSG('k3sqrt','dfunc','dfunc is nan.',-94,0)
endfunction

function d2func(func,r,spiral)
USE STELLARDISK_MODEL
!
! Forward differential
!
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE PRECISION                ::dr
DOUBLE PRECISION,PARAMETER      ::coe(4)=(/2.d0, -5.d0, 4.d0, -1.d0/)
DOUBLE PRECISION                ::d2func,ans,funcs(4)
INTEGER                         ::i
interface 
        function func(x,spiral)
        USE STELLARDISK_MODEL
        DOUBLE PRECISION        ::func
        DOUBLE PRECISION        ::x
        type(typspiral),TARGET                  ::spiral
        ENDFUNCTION func
endinterface        

dr = epsilon(r)**0.3*max(r,epsilon(r))
dr = 1d-6
funcs(1) = func(r,spiral)
funcs(2) = func(r+dr,spiral)
funcs(3) = func(r+2.d0*dr,spiral)
funcs(4) = func(r+3.d0*dr,spiral)
ans = dot_product(funcs,coe)/dr**2
d2func = ans

if(isnan(d2func))CALL XERMSG('k3sqrt','dfunc','dfunc is nan.',-94,0)
endfunction

function cd2func(func,r,spiral)
USE STELLARDISK_MODEL
!
! Forward differential
!
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE PRECISION                ::dr
!DOUBLE PRECISION,PARAMETER      ::coe(4)=(/2.d0, -5.d0, 4.d0, -1.d0/)
DOUBLE COMPLEX,PARAMETER      ::coe(5)=(/35.d0/12.d0,-26.d0/3.d0,19.d0/2.d0,-14.d0/3.d0,11.d0/12.d0/)
DOUBLE COMPLEX                  ::cd2func,ans,funcs(5)
INTEGER                         ::i
interface 
        function func(x,spiral)
        USE STELLARDISK_MODEL
        DOUBLE COMPLEX          ::func
        DOUBLE PRECISION        ::x
        type(typspiral),TARGET                  ::spiral
        ENDFUNCTION func
endinterface        

dr = epsilon(r)**0.5*max(r,epsilon(0.d0))
dr = 1d-3
funcs(1) = func(r,spiral)
funcs(2) = func(r+dr,spiral)
funcs(3) = func(r+2.d0*dr,spiral)
funcs(4) = func(r+3.d0*dr,spiral)
funcs(5) = func(r+4.d0*dr,spiral)
ans = dot_product(dcmplx(coe),funcs)/dr**2
cd2func = ans

if(isnan(real(cd2func)))CALL XERMSG('k3sqrt','dfunc','dfunc is nan.',-94,0)
endfunction
 
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
CALL set_co(spiral)
if(spiral.bndu0)then
        ui = (/dcmplx(a),dcmplx(0.d0),2.d0*sqrt(-q(0.d0))/)
else
        ui = (/dcmplx(a),dcmplx(1.d0,0.d0),dcmplx(0.d0)/)
endif
CALL rk4(a,b,spiral.n,p,q,p,spiral.u,ui)
spiral.ucaled = .true.
spiral.r(:) = real(spiral.u(1,:))
spiral.error = error(spiral) !spiral.fortoone is assigned here
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

FUNCTION CO()
IMPLICIT NONE
DOUBLE PRECISION                        ::CO
DOUBLE PRECISION                        ::B,C,R,RE,AE
INTEGER                                 ::IFLAG

B = 5.d0
C = 14.d0
R = (C-B)/2.d0
RE = 1d-7
AE = 1d-7

CALL DFZERO(fco,B,C,R,RE,AE,IFLAG)
CO = c
ENDFUNCTION

FUNCTION fco(r)
IMPLICIT NONE
DOUBLE PRECISION                ::fco,r
        fco = Omega(r,spiral)-real(spiral.w)/2.d0
ENDFUNCTION

ENDSUBROUTINE

SUBROUTINE set_co(spiral)
type(typspiral),TARGET                  ::spiral

!IF(spiral.co<0.d0)THEN
        spiral.co = CO()
        CALL set_fourtoone(spiral)
!ENDIF
contains

FUNCTION CO()
IMPLICIT NONE
DOUBLE PRECISION                        ::CO
DOUBLE PRECISION                        ::B,C,R,RE,AE
INTEGER                                 ::IFLAG

B = 0.d0
C = 20.d0
R = (C-B)/2.d0
RE = 1d-7
AE = 1d-7

CALL DFZERO(fco,B,C,R,RE,AE,IFLAG)
CO = c
ENDFUNCTION

FUNCTION fco(r)
IMPLICIT NONE
DOUBLE PRECISION                ::fco,r
        fco = Omega(r,spiral)-real(spiral.w)/2.d0
ENDFUNCTION

ENDSUBROUTINE

SUBROUTINE set_fourtoone(spiral)
USE STELLARDISK_MODEL
USE math
IMPLICIT NONE          
type(typspiral)                         ::spiral
DOUBLE COMPLEX          ::error
DOUBLE COMPLEX          ::uu(3)
DOUBLE PRECISION        ::h=10d-5   ,r
DOUBLE PRECISION        ::RE,AE,B,C,RR
INTEGER                 ::l,IFLAG

B = 0.d0
C = 20.d0
RR = 8.d0
RE = 1d-8
AE = 1d-8
call DFZERO(four21,B,C,RR,RE,AE,IFLAG)
r = B
spiral.fortoone = r

CONTAINS

function four21(r)
IMPLICIT NONE
DOUBLE PRECISION                ::four21,r
four21 = (real(spiral.w)- 2.d0*Omega(r,spiral))/kappa(r,spiral) - 0.5d0
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
spiral.h1(1) = (0.d0,0.d0)

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
DOUBLE PRECISION,PARAMETER      ::cutoff = 10.76d0
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

function Sigma(r,spiral) RESULT(ans)
USE STELLARDISK_MODEL
!This is NOT related to density
IMPLICIT NONE
type(typspiral)                         ::spiral
DOUBLE PRECISION                ::ans
DOUBLE PRECISION,INTENT(IN)     ::r
ans = 2.d0*pi*GravConst*sigma0(r,spiral)/snsd(r,spiral)**2
ENDFUNCTION

FUNCTION sigma1(r,th,spiral)
USE STELLARDISK_MODEL
USE math
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE COMPLEX                  ::uu,hh1
DOUBLE PRECISION                ::sigma1
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE PRECISION                ::rad,th
INTEGER                         ::i,j,k,l,n
!interploting u at non-grid point r
n = spiral.n

hh1 = cintplt(spiral.h1,spiral.r,r)

!find density
th = th + spiral.phase
sigma1 = real(hh1*sigma0(r,spiral)/snsd(r,spiral)**2*exp(-2.d0*th*(0.d0,1.d0)))
if(r>spiral.co)sigma1 = sigma1 *exp(-(r-spiral.co)**2/1.0d0**2)

ENDFUNCTION

FUNCTION sigma1r(r,spiral)
USE STELLARDISK_MODEL
USE math
!This is to find density perturbation by solve the k3sqr ODE
IMPLICIT NONE
type(typspiral),TARGET                  ::spiral
DOUBLE PRECISION,INTENT(IN)             ::r
DOUBLE COMPLEX                          ::h1
DOUBLE PRECISION                        ::sigma1r

h1 = cintplt(spiral.h1,spiral.r,r)
sigma1r = abs(h1)/snsd(r,spiral)**2*sigma0(r,spiral)

ENDFUNCTION

function error(spiral)
USE STELLARDISK_MODEL
USE math
IMPLICIT NONE          
type(typspiral)                         ::spiral
DOUBLE COMPLEX          ::error
DOUBLE COMPLEX          ::uu(3)
DOUBLE PRECISION        ::h=2d-1   ,r

r = spiral.fortoone 
uu(2) = cintplt(spiral.u(2,:),spiral.r,r)
uu(3) = cintplt(spiral.u(3,:),spiral.r,r)
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

FUNCTION BulgeSurfaceDensity(r,Spiral)
DOUBLE PRECISION                        ::BulgeSurfaceDensity
DOUBLE PRECISION,INTENT(IN)             ::r
type(typspiral),INTENT(IN),TARGET       ::spiral
DOUBLE PRECISION,POINTER                ::para(:)
DOUBLE PRECISION                        ::BOUND,EPSREL,EPSABS
DOUBLE PRECISION                        ::ans
DOUBLE PRECISION                        ::ABSERR
INTEGER                                 ::NEVAL,IERR,LIMIT,LENW,LAST,INF
DOUBLE PRECISION,ALLOCATABLE            ::WORK(:)
INTEGER,ALLOCATABLE                     ::IWORK(:)
DOUBLE PRECISION                        ::Mb,rb !bulge parameters
DOUBLE PRECISION                        ::rr

para=>spiral.para

rr = r
Mb    = para(3)
rb    = para(4)

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
BulgeSurfaceDensity = ANS
DEALLOCATE(WORK)
DEALLOCATE(IWORK)
CONTAINS 
FUNCTION FUN(z)
DOUBLE PRECISION                        ::FUN,z
        FUN = Mb*(1+(r**2+z**2)/rb**2)**-1.5
ENDFUNCTION
ENDFUNCTION

FUNCTION BulgeSurfaceDensityA(r,Spiral)
DOUBLE PRECISION                        ::BulgeSurfaceDensityA
type(typspiral),INTENT(IN),TARGET       ::spiral
DOUBLE PRECISION,INTENT(IN)             ::r
DOUBLE PRECISION,POINTER                ::para(:)
DOUBLE PRECISION                        ::Mb,rb !bulge parameters
DOUBLE PRECISION                        ::rr
para=>spiral.para

Mb    = para(3)
rb    = para(4)
BulgeSurfaceDensityA = 2.d0*Mb/(1.d0+r**2/rb**2)**1.5*sqrt(r**2+rb**2)
ENDFUNCTION

ENDMODULE STELLARDISK
