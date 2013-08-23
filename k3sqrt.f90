!> @file
!! Main library of StellarSpiral

!> Module for all Omega components.
!> Functions are accessible by function Omega
MODULE ModOmega
USE STELLARDISK_MODEL
IMPLICIT NONE
PRIVATE
PUBLIC                                  ::OmegaOpt
CONTAINS

FUNCTION OmegaOpt(r,spiral,arg)
USE NUM
CLASS(typspiralmodel),INTENT(IN)        ::spiral
REAL(REAL128)                           ::OmegaOpt,Omegasqr
REAL(REAL128),INTENT(IN)                ::r
CHARACTER(*),INTENT(IN),OPTIONAL        ::arg

OmegaOpt = 0
Omegasqr = 0
IF(PRESENT(arg))THEN
        IF(SCAN(arg,'b').ne.0)Omegasqr = Omegasqr + bOmega(r,spiral)**2
        IF(SCAN(arg,'h').ne.0)Omegasqr = Omegasqr + hOmega(r,spiral)**2
        IF(SCAN(arg,'d').ne.0)Omegasqr = Omegasqr + dOmega(r,spiral)**2
        OmegaOpt = sqrt(Omegasqr)
ELSE
        OmegaOpt = TotalOmega(spiral,r)
ENDIF

ENDFUNCTION

FUNCTION TotalOmega(spiral,r)
USE NUM
IMPLICIT NONE
REAL(REAL128)                           ::TotalOmega
REAL(REAL128),INTENT(IN)                ::r
CLASS(typspiralmodel)                   ::spiral

TotalOmega = bOmega(r,spiral) &
           + dOmega(r,spiral) &
           + hOmega(r,spiral)

ENDFUNCTION

FUNCTION bOmega(r,spiral)
USE NUM
IMPLICIT NONE
REAL(REAL128)          ::bOmega
REAL(REAL128),INTENT(IN)::r
!bulge
REAL(REAL128)          ::rb,Mb,gBulge,VBulge
!spiral
type(typspiralmodel),TARGET                  ::spiral
REAL(REAL128),POINTER                ::para(:)=>null()

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
gBulge = (-r/sqrt(1.d0+r**2/rb**2)/rb+asinh(r/rb))/r**2

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
REAL(REAL128)          ::dOmega
REAL(REAL128),INTENT(IN)::r
!disk
REAL(REAL128)          ::dM,da,db
DOUBLE COMPLEX            ::VDisk
REAL(REAL128)          ::a1,a2,M1,M2
!spiral
type(typspiralmodel),TARGET                  ::spiral
REAL(REAL128),POINTER                ::para(:)=>null()

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
REAL(REAL128)        ::pDisk,r
pDisk  = -GravConst*dM
pDisk  = pDisk/sqrt(r**2+(da+db)**2)
ENDFUNCTION

FUNCTION dpDisk(r)
IMPLICIT NONE
REAL(REAL128)        ::dpDisk
REAL(REAL128)        ::r
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
REAL(REAL128)                ::r
REAL(REAL128)                ::x1,x2
REAL(REAL128)                ::a1,a2,M1,M2

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
REAL(REAL128)                ::H,x

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
type(typspiralmodel),TARGET             ::spiral
REAL(REAL128)                           ::hOmega
REAL(REAL128),INTENT(IN)                ::r
!Halo
REAL(REAL128)                           ::Lh,rhoh,gHalo,VHalo
!spiral
REAL(REAL128),POINTER                   ::para(:)=>null()

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
ENDMODULE

MODULE STELLARDISK
USE STELLARDISK_MODEL
USE ModOmega
IMPLICIT NONE
TYPE,EXTENDS(typspiralmodel)            ::typspiral
        CONTAINS
        PROCEDURE                       ::Omega
ENDTYPE

CONTAINS

FUNCTION typspiral_new(n,rmax,para,mode)
TYPE(typspiral)                         ::typspiral_new
TYPE(typgalaxy_para)                    ::para
REAL(REAL128)                           ::rmax
INTEGER                                 ::n,mode

CALL typspiral_new.init(n,rmax,para,mode)

ENDFUNCTION

!> Omega
!> @param[in] r REAL
!! @param[in] arg CHARACTER, Requested component(s).
! @param[out] REAL128
FUNCTION Omega(spiral,r,arg)
IMPLICIT NONE
CLASS(typspiral)                        ::spiral
REAL(REAL128)                           ::Omega
REAL(REAL128)                           ::r
CHARACTER(*),OPTIONAL                   ::arg
IF(PRESENT(arg))THEN
        Omega = OmegaOpt(r,spiral,arg)
ELSE
        Omega = OmegaOpt(r,spiral)
ENDIF
ENDFUNCTION

ENDMODULE STELLARDISK
