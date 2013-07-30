subroutine fluxHLLAdiMHD(ql,qr,slopeL,slopeR,Ftemp)
use common_params
implicit none
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),Ftemp(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spd ! sound speed
double precision:: rhoL,vxL,vyL,vzL,bxL,byL,bzL,eneL,sndfL,pL,BsqL,vsqL,pTL
double precision:: rhoR,vxR,vyR,vzR,bxR,byR,bzR,eneR,sndfR,pR,BsqR,vsqR,pTR
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
vzL =UL(4)
bxL =0.5d0*(UL(5)+UR(5))
byL =UL(6)
bzL =UL(7)
pL = UL(8)
!eneL = UL(8)
vsqL = vxL**2.d0+vyL**2.d0+vzL**2.d0
BsqL = bxL**2.d0+byL**2.d0+bzL**2.d0
eneL = (pL/(gam-1.d0)+0.5d0*rhoL*vsqL+0.5d0*BsqL)/rhoL
!pL = (gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*vsqL-0.5d0*BsqL) ! gaseous pressure
sndfL = dsqrt((gam*pL+BsqL+dsqrt((gam*pL+BsqL)**2.d0-4.d0*gam*pL*bxL**2.d0))/(2.d0*rhoL)) ! fast wave
pTL = pL+0.5d0*BsqL

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
vzR =UR(4)
bxR =0.5d0*(UL(5)+UR(5))
byR =UR(6)
bzR =UR(7)
!eneR=UR(8)
pR=UR(8)
vsqR = vxR**2.d0+vyR**2.d0+vzR**2.d0
BsqR = bxR**2.d0+byR**2.d0+bzR**2.d0
eneR = (pR/(gam-1.d0)+0.5d0*rhoR*vsqR+0.5d0*BsqR)/rhoR
!pR = (gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*vsqR-0.5d0*BsqR) ! gaseous pressure
sndfR = dsqrt((gam*pR+BsqR+dsqrt((gam*pR+BsqR)**2.d0-4.d0*gam*pR*bxR**2.d0))/(2.d0*rhoR)) ! fast wave
pTR = pR+0.5d0*BsqR

SL=dmin1(vxL,vxR)-dmax1(sndfL,sndfR)
SR=dmax1(vxL,vxR)+dmax1(sndfL,sndfR)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL-bxL**2.d0
FL(3)=rhoL*vxL*vyL-bxL*byL
FL(4)=rhoL*vxL*vzL-bxL*bzL
FL(5)=0.d0
FL(6)=byL*vxL-bxL*vyL
FL(7)=bzL*vxL-bxL*vzL
FL(8)=(rhoL*eneL+pTL)*vxL-bxL*(vxL*bxL+vyL*byL+vzL*bzL)

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR-bxR**2.d0
FR(3)=rhoR*vxR*vyR-bxR*byR
FR(4)=rhoR*vxR*vzR-bxR*bzR
FR(5)=0.d0
FR(6)=byR*vxR-bxR*vyR
FR(7)=bzR*vxR-bxR*vzR
FR(8)=(rhoR*eneR+pTR)*vxR-bxR*(vxR*bxR+vyR*byR+vzR*bzR)

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)
Fs(4)=(SR*FL(4)-SL*FR(4)+SR*SL*(UR(1)*UR(4)-UL(1)*UL(4)))/(SR-SL)
Fs(5)=0.d0
Fs(6)=(SR*FL(6)-SL*FR(6)+SR*SL*(UR(6)-UL(6)))/(SR-SL)
Fs(7)=(SR*FL(7)-SL*FR(7)+SR*SL*(UR(7)-UL(7)))/(SR-SL)
Fs(8)=(SR*FL(8)-SL*FR(8)+SR*SL*(UR(1)*eneR-UL(1)*eneL))/(SR-SL)

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
