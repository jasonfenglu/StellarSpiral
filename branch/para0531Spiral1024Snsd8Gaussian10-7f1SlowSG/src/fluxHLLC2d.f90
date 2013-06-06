subroutine fluxHLLC2d(ql,qr,slopeL,slopeR,Ftemp)
use common_params
implicit none
double precision::Ftemp(nvar),ql(nvar),qr(nvar),slopeL(nvar),slopeR(nvar), &
                  FL(nvar),FR(nvar),FsL(nvar), FsR(nvar),Fs(nvar)
double precision::rhoL,rhoR,vxL,vxR,seneL,seneR,eneL,eneR,pTL,pTR,pxL,pxR, &
             spdL,spdR,SL,SR,SM,pstar,rhosL,rhosR,pxsL,pxsR,enesL,enesR, &
             vyL,vyR,pyL,pyR,pysL,pysR,pL,pR

rhoL =ql(1)+0.5d0*slopeL(1)
vxL  =ql(2)+0.5d0*slopeL(2)
vyL  =ql(3)+0.5d0*slopeL(3)
pL   =ql(4)+0.5d0*slopeL(4)
eneL =pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2)
!seneL=ql(4)+0.5d0*slopeL(4)
!eneL =seneL*rhoL

rhoR =qr(1)-0.5d0*slopeR(1)
vxR  =qr(2)-0.5d0*slopeR(2)
vyR  =qr(3)-0.5d0*slopeR(3)
pR   =qr(4)-0.5d0*slopeR(4)
eneR =pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2)
!seneR=qr(4)-0.5d0*slopeR(4)
!eneR =seneR*rhoR

pTL = (gam-1.d0)*(eneL-0.5d0*rhoL*(vxL**2+vyL**2))
pTR = (gam-1.d0)*(eneR-0.5d0*rhoR*(vxR**2+vyR**2))

!if(pTL .le. 0.d0 .or. pTR .le. 0.d0) then
if(pL .le. 0.d0 .or. pR .le. 0.d0) then
rhoL =ql(1)
vxL  =ql(2)
vyL  =ql(3) 
pL   =ql(4)
eneL =pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2)
!seneL=ql(4)
!eneL =seneL*rhoL

rhoR =qr(1)
vxR  =qr(2)
vyR  =qr(3)
pR   =qr(4)
eneR =pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2)
!seneR=qr(4)
!eneR =seneR*rhoR

pTL = (gam-1.d0)*(eneL-0.5d0*rhoL*(vxL**2+vyL**2))
pTR = (gam-1.d0)*(eneR-0.5d0*rhoR*(vxR**2+vyR**2))
!write(*,*) "negative pressure!"
endif

spdL = dsqrt(gam*pTL/rhoL)
pxL = vxL*rhoL
pyL = vyL*rhoL

spdR = dsqrt(gam*pTR/rhoR)
pxR = vxR*rhoR
pyR = vyR*rhoR

SL = dmin1(vxL-spdL,vxR-spdR)
SR = dmax1(vxL+spdL,vxR+spdR)
SM = ((SR-vxR)*rhoR*vxR-(SL-vxL)*rhoL*vxL-pTR+pTL)  &
        /((SR-vxR)*rhoR-(SL-vxL)*rhoL)
pstar = rhoL*(vxL-SL)*(vxL-SM)+pTL

rhosL = rhoL*(SL-vxL)/(SL-SM)
pxsL = ((SL-vxL)*pxL+(pstar-pTL))/(SL-SM)
pysL = rhosL*vyL
enesL = ((SL-vxL)*eneL-pTL*vxL+pstar*SM)/(SL-SM)

rhosR = rhoR*(SR-vxR)/(SR-SM)
pxsR = ((SR-vxR)*pxR+(pstar-pTR))/(SR-SM)
pysR = rhosR*vyR
enesR = ((SR-vxR)*eneR-pTR*vxR+pstar*SM)/(SR-SM)

FL(1) = rhoL*vxL
FL(2) = rhoL*vxL**2.d0 +pTL
FL(3) = rhoL*vxL*vyL
FL(4) = (eneL+pTL)*vxL

FR(1) = rhoR*vxR
FR(2) = rhoR*vxR**2.d0 +pTR
FR(3) = rhoR*vxR*vyR
FR(4) = (eneR+pTR)*vxR

FsL(1)= FL(1) + SL*(rhosL-rhoL)
FsL(2)= FL(2) + SL*(pxsL-pxL)
FsL(3)= FL(3) + SL*(pysL-pyL)
FsL(4)= FL(4) + SL*(enesL-eneL)

FsR(1)= FR(1) + SR*(rhosR-rhoR)
FsR(2)= FR(2) + SR*(pxsR-pxR)
FsR(3)= FR(3) + SR*(pysR-pyR)
FsR(4)= FR(4) + SR*(enesR-eneR)

IF (SL .gt. 0.d0)then
       Ftemp=FL
ELSEIF (SL .le. 0.d0 .and. SM .ge. 0.d0)then
       Ftemp=FsL
ELSEIF (SM .le. 0.d0 .and. SR .ge. 0.d0)then
       Ftemp=FsR
ELSEIF (SR .lt. 0.d0) then
       Ftemp=FR
ENDIF

end subroutine
