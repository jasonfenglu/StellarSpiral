subroutine fluxRieIso3d(ql,qr,slopeL,slopeR,Ftemp)
use common_params
implicit none
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),Ftemp(NVAR),qm(NVAR)
double precision:: rhoL,rhoR,vxL,vxR,vyL,vyR,vzL,vzR
double precision::c,shkR,rarR,csq,vstar,rhostar,rho0,rhosR,srhosR,aisrhosR,sigma2, &
                   rhosL,srhosL,aisrhosL,sigma1,aa,bb,uu,srhostar
integer:: Intr,n
double precision:: tolerr = 1.d-10

ql(1) = ql(1)+0.5d0*slopeL(1)
ql(2) = ql(2)+0.5d0*slopeL(2)
ql(3) = ql(3)+0.5d0*slopeL(3)
ql(4) = ql(4)+0.5d0*slopeL(4)

qr(1) = qr(1)-0.5d0*slopeR(1)
qr(2) = qr(2)-0.5d0*slopeR(2)
qr(3) = qr(3)-0.5d0*slopeR(3)
qr(4) = qr(4)-0.5d0*slopeR(4)

c = snd
Intr = 10


if (ql(1).eq.qr(1).and.ql(2).eq.qr(2)) then
     qm(1) = qr(1)
     qm(2) = qr(2)
if(ql(2).ge. 0.   ) then
     qm(3) = ql(3)
     qm(4) = ql(4)
else
     qm(3) = qr(3)
     qm(4) = qr(4)
endif
goto 1000
endif

      rhoL = ql(1)
      rhoR = qr(1)
       vxL = ql(2)
       vxR = qr(2)
       vyL = ql(3)
       vyR = qr(3)
       vzL = ql(4)
       vzR = qr(4)


      shkR = -c*(rhoR-rhoL)/dsqrt(rhoR*rhoL)
      rarR = -c*dlog(rhoR/rhoL)
       csq =  c*c

      if(rhoL .ge. rhoR) then
        if ((vxR-vxL) .gt. rarR) then
          goto 100
        elseif ((vxR-vxL) .gt. -shkR) then
          goto 200
        else
          goto 400
        endif
      else 
        if ((vxR-vxL) .lt. shkR ) then
          goto 400
        elseif ((vxR-vxL) .lt. -rarR) then
          goto 300
        else
          goto 100
        endif
      endif

!cccccc 1-rarefaction wave, 2-rarefaction wave
!cccccc 1-rarefaction wave, 2-rarefaction wave
!cccccc 1-rarefaction wave, 2-rarefaction wave
100   vstar   = 0.5d0*(vxL+vxR)+0.5d0*c*dlog(rhoL/rhoR)
      rhostar = rhoL*dexp(-(vstar-vxL)/c)

      if ((vxL-c) .ge. 0.) then
         qm(2) = vxL
         qm(1) = rhoL
      elseif ((vstar-c) .ge. 0.) then
         qm(2) =  c
         qm(1) =  rhoL*dexp((vxL-c)/c)
      elseif ((vstar+c) .ge. 0.) then
         qm(2) =  vstar
         qm(1) =  rhostar
      elseif ((  vxR+c) .ge. 0.) then
         qm(2) =  -c
         qm(1) =  rhoR*dexp(-(vxR+c)/c)
      else
         qm(2) =   vxR
         qm(1) =  rhoR
      endif
      if ( qm(2) .ge. 0.) then
         qm(3) =  ql(3)
         qm(4) =  ql(4)
      else
         qm(3) =  qr(3)
         qm(4) =  qr(4)
      endif
      goto 1000

!cccccc1-rarefaction wave, 2-shock wave
!cccccc1-rarefaction wave, 2-shock wave
!cccccc1-rarefaction wave, 2-shock wave
200   rho0 = rhoR
      do n = 1, INtr
        rhosR    = rho0/rhoR
        srhosR   = dsqrt(rhosR)
        aisrhosR = 1.d0/srhosR
        rhostar  = rho0 - &
       (vxR-vxL+c*(srhosR-aisrhosR)+c*dlog(rho0/rhoL)) &
       /(0.5d0*c/rho0*(srhosR+aisrhosR)+c/rho0)
        if (dabs((rho0-rhostar)/rho0) .lt. tolerr) then
            goto 210
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1R2S'
      write(6,*) 'rhostar,rhoL=',rhostar,rhoL
      write(6,*) 'uxL    , uxR=',vxL,      vxR
      stop
210   vstar  = vxL - c*dlog(rhostar/rhoL)
      sigma2 = vxR + c*dsqrt(rhostar/rhoR)
      if ( sigma2 .le. 0.) then
         qm(2) =  vxR
         qm(1) = rhoR
      elseif ((vstar-c) .le. 0.) then
         qm(2) = vstar
         qm(1) = rhostar
      elseif ((vxL-c) .le. 0.) then
         qm(2) =  c
         qm(1) =  rhoL*dexp((vxL-c)/c)
      else
         qm(2) = vxL
         qm(1) = rhoL
      endif
      if ( qm(2) .ge. 0.) then
         qm(3) = ql(3)
         qm(4) = ql(4)
      else
         qm(3) = qr(3)
         qm(4) = qr(4)
      endif
      goto 1000
!cccccc 1-shock wave, 2-rarefaction wave
!cccccc 1-shock wave, 2-rarefaction wave
!cccccc 1-shock wave, 2-rarefaction wave
300   rho0 = rhoL
      do n = 1, INtr
        rhosL    = rho0/rhoL
        srhosL   = dsqrt(rhosL)
        aisrhosL = 1.d0/srhosL
        rhostar  = rho0- &
      (vxR-vxL+c*(srhosL-aisrhosL)+c*dlog(rho0/rhoR)) &
       /(0.5d0*c/rho0*(srhosL+aisrhosL)+c/rho0)
        if (dabs((rho0-rhostar)/rho0) .lt. tolerr) then
            goto 310
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1S2R'
      write(6,*) 'rhostar,rhoL=',rhostar,rhoL
      write(6,*) 'uxL    , uxR=',vxL,      vxR
      stop
310   vstar  = vxR + c*dlog(rhostar/rhoR)
      sigma1 = vxL - c*dsqrt(rhostar/rhoL)
      if (sigma1 .ge. 0.) then
         qm(2) = vxL
         qm(1) = rhoL
      elseif ((vstar+c) .ge. 0.) then
         qm(2) = vstar
         qm(1) = rhostar
      elseif ((vxR+c) .ge. 0.) then
         qm(2) =  -c
         qm(1) =  rhoR*dexp(-(vxR+c)/c)
      else
         qm(2) = vxR
         qm(1) = rhoR
      endif
      if ( qm(2) .ge. 0.) then
         qm(3) = ql(3)
         qm(4) = ql(4)
      else
         qm(3) = qr(3)
         qm(4) = qr(4)
      endif
      goto 1000
!       1-shock wave, 2-shock wave
!       1-shock wave, 2-shock wave
!       1-shock wave, 2-shock wave
400   aa       = 1.d0/dsqrt(rhoL)+1.d0/dsqrt(rhoR)
      bb       =    dsqrt(rhoL)+   dsqrt(rhoR)
      uu       = vxR-vxL
      srhostar = (-uu+dsqrt(uu**2+4.d0*csq*aa*bb))/(2.d0*c*aa)
      rhostar  = srhostar*srhostar
      vstar    = vxL - c*(dsqrt(rhostar/rhoL)-dsqrt(rhoL/rhostar))
      sigma1   = vxL - c*dsqrt(rhostar/rhoL)
      sigma2   = vxR + c*dsqrt(rhostar/rhoR)
      if     (sigma1 .ge. 0.) then
       qm(2) = vxL
       qm(1) = rhoL
      elseif (sigma2 .ge. 0.) then
       qm(2) = vstar
       qm(1) = rhostar
      else
       qm(2) = vxR
       qm(1) = rhoR
      endif
         if  (qm(2) .ge. 0.) then
       qm(3) = ql(3)
       qm(4) = ql(4)
      else
       qm(3) = qr(3)
       qm(4) = qr(4)
      endif
      goto 1000

1000  Ftemp(1) = qm(1)*(   qm(2)          )
      Ftemp(2) = qm(1)*( qm(2)*qm(2) + c*c)
      Ftemp(3) = qm(1)*( qm(2)*qm(3)      )
      Ftemp(4) = qm(1)*( qm(2)*qm(4)      )
      Ftemp(5) = 0.d0
end subroutine
