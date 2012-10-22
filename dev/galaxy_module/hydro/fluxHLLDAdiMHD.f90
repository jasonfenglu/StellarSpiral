subroutine fluxHLLDAdiMHD(ql,qr,slopeL,slopeR,Ftemp)
use common_params
implicit none
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),Ftemp(NVAR)
double precision::rhoL,vxL,vyL,vzL,bxL,byL,bzL,eneL,vsqL,BsqL,pL,pTL,fastL
double precision::rhoR,vxR,vyR,vzR,bxR,byR,bzR,eneR,vsqR,BsqR,pR,pTR,fastR
double precision::SM,SL,SR,SLs,SRs,pTs
double precision::vxLs,vxLss,vxRs,vxRss,vyss,vzss
double precision::pTLs,pTLss,pTRs,pTRss,byss,bzss
double precision::rhoLs,rhoRs

double precision,dimension(:),allocatable::UL,UR,ULs,URs,ULss,URss
double precision,dimension(:),allocatable::FL,FR,FLs,FRs,FLss,FRss
double precision:: tempsL,tempsR,tempss

allocate(UL(NVAR))
allocate(UR(NVAR))
allocate(ULs(NVAR))
allocate(URs(NVAR))
allocate(ULss(NVAR))
allocate(URss(NVAR))
allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(FLs(NVAR))
allocate(FRs(NVAR))
allocate(FLss(NVAR))
allocate(FRss(NVAR))

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
!eneL=UL(8)
vsqL=vxL**2.d0+vyL**2.d0+vzL**2.d0
BsqL=bxL**2.d0+byL**2.d0+bzL**2.d0
eneL = (pL/(gam-1.d0)+0.5d0*rhoL*vsqL+0.5d0*BsqL)/rhoL
!pL  =(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*vsqL-0.5d0*BsqL)
pTL =pL+0.5d0*BsqL
fastL=dsqrt((gam*pL+BsqL+dsqrt((gam*pL+BsqL)**2.d0-4.d0*gam*pL*bxL**2))/(2.d0*rhoL))

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
vzR =UR(4)
bxR =0.5d0*(UL(5)+UR(5))
byR =UR(6)
bzR =UR(7)
!eneR=UR(8)
pR = UR(8)
vsqR=vxR**2.d0+vyR**2.d0+vzR**2.d0
BsqR=bxR**2.d0+byR**2.d0+bzR**2.d0
eneR = (pR/(gam-1.d0)+0.5d0*rhoR*vsqR+0.5d0*BsqR)/rhoR
!pR  =(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*vsqR-0.5d0*BsqR)
pTR =pR+0.5d0*BsqR
fastR=dsqrt((gam*pR+BsqR+dsqrt((gam*pR+BsqR)**2.d0-4.d0*gam*pR*bxR**2))/(2.d0*rhoR))

SL=dmin1(vxL,vxR)-dmax1(fastL,fastR)
SR=dmax1(vxL,vxR)+dmax1(fastL,fastR)

SM=((SR-vxR)*rhoR*vxR-(SL-vxL)*rhoL*vxL-pTR+pTL)/((SR-vxR)*rhoR-(SL-vxL)*rhoL)
pTs=((SR-vxR)*rhoR*pTL-(SL-vxL)*rhoL*pTR+rhoL*rhoR*(SR-vxR)*(SL-vxL)*(vxR-vxL))/((SR-vxR)*rhoR-(SL-vxL)*rhoL)

vxLs =SM
vxLss=SM
vxRss=SM
vxRs =SM

pTLs =pTs
pTLss=pTs
pTRs =pTs
pTRss=pTs

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL-bxL**2.d0
FL(3)=rhoL*vyL*vxL-bxL*byL
FL(4)=rhoL*vzL*vxL-bxL*bzL
FL(5)=0.d0
FL(6)=byL*vxL-bxL*vyL
FL(7)=bzL*vxL-bxL*vzL
FL(8)=(rhoL*eneL+pTL)*vxL-bxL*(vxL*bxL+vyL*byL+vzL*bzL)

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR-bxR**2.d0
FR(3)=rhoR*vyR*vxR-bxR*byR
FR(4)=rhoR*vzR*vxR-bxR*bzR
FR(5)=0.d0
FR(6)=byR*vxR-bxR*vyR
FR(7)=bzR*vxR-bxR*vzR
FR(8)=(rhoR*eneR+pTR)*vxR-bxR*(vxR*bxR+vyR*byR+vzR*bzR)

! ============================
if(SL .eq. SM) then
ULs(1)=rhoL
else
ULs(1)=rhoL*(SL-vxL)/(SL-SM)
endif
ULs(2)=vxLs
tempsL=rhoL*(SL-vxL)*(SL-SM)-bxL**2.d0
if(tempsL .eq. 0.d0) then
ULs(3)=vyL
ULs(4)=vzL
ULs(5)=bxL
ULs(6)=byL
ULs(7)=bzL
else
ULs(3)=vyL-bxL*byL*(SM-vxL)/tempsL
ULs(4)=vzL-bxL*bzL*(SM-vxL)/tempsL
ULs(5)=bxL
ULs(6)=byL*(rhoL*(SL-vxL)**2.d0-bxL**2)/tempsL
ULs(7)=bzL*(rhoL*(SL-vxL)**2.d0-bxL**2)/tempsL
endif
ULs(8)=((SL-vxL)*rhoL*eneL-pTL*vxL+pTs*SM+bxL*(vxL*bxL+vyL*byL+vzL*bzL-ULs(2)*ULs(5)-ULs(3)*ULs(6)-ULs(4)*ULs(7)))/(SL-SM)/ULs(1)

if(SR .eq. SM) then
URs(1)=rhoR
else
URs(1)=rhoR*(SR-vxR)/(SR-SM)
endif
URs(2)=vxRs
tempsR=rhoR*(SR-vxR)*(SR-SM)-bxR**2.d0
if(tempsL .eq. 0.d0) then
URs(3)=vyR
URs(4)=vzR
URs(5)=bxR
URs(6)=byR
URs(7)=bzR
else
URs(3)=vyR-bxR*byR*(SM-vxR)/tempsR
URs(4)=vzR-bxR*bzR*(SM-vxR)/tempsR
URs(5)=bxR
URs(6)=byR*(rhoR*(SR-vxR)**2.d0-bxR**2)/tempsR
URs(7)=bzR*(rhoR*(SR-vxR)**2.d0-bxR**2)/tempsR
endif
URs(8)=((SR-vxR)*rhoR*eneR-pTR*vxR+pTs*SM+bxR*(vxR*bxR+vyR*byR+vzR*bzR-URs(2)*URs(5)-URs(3)*URs(6)-URs(4)*URs(7)))/(SR-SM)/URs(1)
!==========================
tempss=dsqrt(ULs(1))+dsqrt(URs(1))
vyss=(dsqrt(ULs(1))*ULs(3)+dsqrt(URs(1))*URs(3)+(URs(6)-ULs(6))*dsign(1.d0,bxL))/tempss
vzss=(dsqrt(ULs(1))*ULs(4)+dsqrt(URs(1))*URs(4)+(URs(7)-ULs(7))*dsign(1.d0,bxL))/tempss
byss=(dsqrt(ULs(1))*URs(6)+dsqrt(URs(1))*ULs(6)+dsqrt(URs(1)*ULs(1))*(URs(3)-ULs(3))*dsign(1.d0,bxL))/tempss
bzss=(dsqrt(ULs(1))*URs(7)+dsqrt(URs(1))*ULs(7)+dsqrt(URs(1)*ULs(1))*(URs(4)-ULs(4))*dsign(1.d0,bxL))/tempss

ULss(1)=ULs(1)
ULss(2)=vxLss
if(BxL .eq. 0.d0) then
ULss(3)=ULs(3)
ULss(4)=ULs(4)
ULss(5)=0.d0
ULss(6)=ULs(6)
ULss(7)=ULs(7)
ULss(8)=ULs(8)
else
ULss(3)=vyss
ULss(4)=vzss
ULss(5)=bxL
ULss(6)=byss
ULss(7)=bzss
ULss(8)=(ULs(8)*ULs(1)-dsqrt(ULs(1))*(ULs(2)*ULs(5)+ULs(3)*ULs(6)+ULs(4)*ULs(7)-ULss(2)*ULss(5)-ULss(3)*ULss(6)-ULss(4)*ULss(7))*dsign(1.d0,bxL))/ULss(1)
endif

URss(1)=URs(1)
URss(2)=vxRss
if(BxR .eq. 0.d0) then
URss(3)=URs(3)
URss(4)=URs(4)
URss(5)=0.d0
URss(6)=URs(6)
URss(7)=URs(7)
URss(8)=URs(8)
else
URss(3)=vyss
URss(4)=vzss
URss(5)=bxR
URss(6)=byss
URss(7)=bzss
URss(8)=(URs(8)*URs(1)+dsqrt(URs(1))*(URs(2)*URs(5)+URs(3)*URs(6)+URs(4)*URs(7)-URss(2)*URss(5)-URss(3)*URss(6)-URss(4)*URss(7))*dsign(1.d0,bxL))/URss(1)
endif

FLs(1)=FL(1)+SL*(ULs(1)-UL(1))
FLs(2)=FL(2)+SL*(ULs(1)*ULs(2)-UL(1)*UL(2))
FLs(3)=FL(3)+SL*(ULs(1)*ULs(3)-UL(1)*UL(3))
FLs(4)=FL(4)+SL*(ULs(1)*ULs(4)-UL(1)*UL(4))
FLs(5)=FL(5)+SL*(ULs(5)-UL(5))
FLs(6)=FL(6)+SL*(ULs(6)-UL(6))
FLs(7)=FL(7)+SL*(ULs(7)-UL(7))
FLs(8)=FL(8)+SL*(ULs(1)*ULs(8)-UL(1)*eneL)

FRs(1)=FR(1)+SR*(URs(1)-UR(1))
FRs(2)=FR(2)+SR*(URs(1)*URs(2)-UR(1)*UR(2))
FRs(3)=FR(3)+SR*(URs(1)*URs(3)-UR(1)*UR(3))
FRs(4)=FR(4)+SR*(URs(1)*URs(4)-UR(1)*UR(4))
FRs(5)=FR(5)+SR*(URs(5)-UR(5))
FRs(6)=FR(6)+SR*(URs(6)-UR(6))
FRs(7)=FR(7)+SR*(URs(7)-UR(7))
FRs(8)=FR(8)+SR*(URs(1)*URs(8)-UR(1)*eneR)

SLs=SM-dabs(bxL)/dsqrt(ULs(1))
SRs=SM+dabs(bxR)/dsqrt(URs(1))

FLss(1)=FL(1)+SLs*ULss(1)-(SLs-SL)*ULs(1)-SL*UL(1)
FLss(2)=FL(2)+SLs*(ULss(1)*ULss(2))-(SLs-SL)*(ULs(1)*ULs(2))-SL*(UL(1)*UL(2))
FLss(3)=FL(3)+SLs*(ULss(1)*ULss(3))-(SLs-SL)*(ULs(1)*ULs(3))-SL*(UL(1)*UL(3))
FLss(4)=FL(4)+SLs*(ULss(1)*ULss(4))-(SLs-SL)*(ULs(1)*ULs(4))-SL*(UL(1)*UL(4))
FLss(5)=FL(5)+SLs*ULss(5)-(SLs-SL)*ULs(5)-SL*UL(5)
FLss(6)=FL(6)+SLs*ULss(6)-(SLs-SL)*ULs(6)-SL*UL(6)
FLss(7)=FL(7)+SLs*ULss(7)-(SLs-SL)*ULs(7)-SL*UL(7)
FLss(8)=FL(8)+SLs*(ULss(1)*ULss(8))-(SLs-SL)*(ULs(1)*ULs(8))-SL*(UL(1)*eneL)

FRss(1)=FR(1)+SRs*URss(1)-(SRs-SR)*URs(1)-SR*UR(1)
FRss(2)=FR(2)+SRs*(URss(1)*URss(2))-(SRs-SR)*(URs(1)*URs(2))-SR*(UR(1)*UR(2))
FRss(3)=FR(3)+SRs*(URss(1)*URss(3))-(SRs-SR)*(URs(1)*URs(3))-SR*(UR(1)*UR(3))
FRss(4)=FR(4)+SRs*(URss(1)*URss(4))-(SRs-SR)*(URs(1)*URs(4))-SR*(UR(1)*UR(4))
FRss(5)=FR(5)+SRs*URss(5)-(SRs-SR)*URs(5)-SR*UR(5)
FRss(6)=FR(6)+SRs*URss(6)-(SRs-SR)*URs(6)-SR*UR(6)
FRss(7)=FR(7)+SRs*URss(7)-(SRs-SR)*URs(7)-SR*UR(7)
FRss(8)=FR(8)+SRs*(URss(1)*URss(8))-(SRs-SR)*(URs(1)*URs(8))-SR*(UR(1)*eneR)

if(SL .gt. 0.d0) then
  Ftemp=FL
elseif (SL .le. 0.d0 .and. SLs .ge. 0.d0) then
  Ftemp=FLs
elseif (SLs .le. 0.d0 .and. SM .ge. 0.d0) then
  Ftemp=FLss
elseif (SM .le. 0.d0 .and. SRs .ge. 0.d0) then
  Ftemp=FRss
elseif (SRs .le. 0.d0 .and. SR .ge. 0.d0) then
  Ftemp=FRs
elseif (SR .lt. 0.d0) then
  Ftemp=FR
endif

deallocate(UL)
deallocate(UR)
deallocate(ULs)
deallocate(URs)
deallocate(ULss)
deallocate(URss)
deallocate(FL)
deallocate(FR)
deallocate(FLs)
deallocate(FRs)
deallocate(FLss)
deallocate(FRss)
end subroutine
