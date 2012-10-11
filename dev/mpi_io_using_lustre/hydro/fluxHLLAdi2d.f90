subroutine fluxHLLAdi2d(ql,qr,slopeL,slopeR,Ftemp)
use common_params
implicit none
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),Ftemp(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spd ! sound speed
double precision:: rhoL,vxL,vyL,eneL,spdL,pTL,pL
double precision:: rhoR,vxR,vyR,eneR,spdR,pTR,pR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
pL  =UL(4) 
eneL = (pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2))/rhoL
pTL = rhoL*(eneL-0.5d0*(vxL**2.d0+vyL**2.d0))*(gam-1.d0)
spdL = dsqrt(gam*pTL/rhoL)

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
pR  =UR(4)
eneR = (pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2))/rhoR
pTR = rhoR*(eneR-0.5d0*(vxR**2.d0+vyR**2.d0))*(gam-1.d0)
spdR = dsqrt(gam*pTR/rhoR)

SL=dmin1(vxL-spdL,vxR-spdR)
SR=dmax1(vxL+spdL,vxR+spdR)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL
FL(3)=rhoL*vxL*vyL
FL(4)=(rhoL*eneL+pTL)*vxL

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR
FR(3)=rhoR*vxR*vyR
FR(4)=(rhoR*eneR+pTR)*vxR

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)
Fs(4)=(SR*FL(4)-SL*FR(4)+SR*SL*(UR(1)*UR(4)-UL(1)*UL(4)))/(SR-SL)

if(SL .gt. 0.d0) then
   Ftemp=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   Ftemp=Fs
else
   Ftemp=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)
end subroutine
