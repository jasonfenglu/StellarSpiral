subroutine fluxHLLIso3d(ql,qr,slopeL,slopeR,Ftemp)
use common_params
implicit none
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),Ftemp(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spd ! sound speed
double precision:: rhoL,vxL,vyL,vzL,pTL
double precision:: rhoR,vxR,vyR,vzR,pTR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

spd=snd
rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
vzL =UL(4)
pTL =rhoL*spd**2.d0

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
vzR =UR(4)
pTR =rhoR*spd**2.d0

SL=dmin1(vxL-spd,vxR-spd)
SR=dmax1(vxL+spd,vxR+spd)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL
FL(3)=rhoL*vxL*vyL
FL(4)=rhoL*vxL*vzL
FL(5)=0.d0

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR
FR(3)=rhoR*vxR*vyR
FR(4)=rhoR*vxR*vzR
FR(5)=0.d0

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)
Fs(4)=(SR*FL(4)-SL*FR(4)+SR*SL*(UR(1)*UR(4)-UL(1)*UL(4)))/(SR-SL)
Fs(5)=0.d0

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
