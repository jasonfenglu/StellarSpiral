c THIS PROGRAM IS GENERATED AUTOMATICALLY BY LAMAC PROGRAM GENERATOR
C WEB INFORMATION http://140.136.192.14/~yen/lamac/lamac2d.php
C Version 1.0
C Please do not modify it.
C If you have any problem, please contact with
C lamac@asiaa.sinica.edu.tw
C Date of Generateion :7/21/2004 Time:14:5:34(Hours:Minutes:Seconds)
C Informations of Program
C Project Name = spiral
C23456789212345678903123456789412345678951234567896123456789712
C program (begin)
      USE KEPLER
      USE find_dt
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
ccccc MPI parallel computation
      include 'mpif.h'
ccccc set global parameters
      include 'parameters'
ccccc q,qa are (rho,rho*ux,rho*uy)=(density,x-momentum,y-momentum)
      dimension  q(ibeg-ibuf  :iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
      dimension qa(ibeg-ibuf  :iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
ccccc MPI parallel computation initialization -start*
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      write(6,*) nprocs, npc, myid
      if (nprocs.ne.npc) call MPI_Abort(ierr)
ccccc MPI parallel computation initialization -stop*
ccccc get variables
      call vars(xmin,xmax,ymin,ymax,
     &          dcfl,ncir,psrl,ichk,
     &          fstr,fstp,nstp,unk1)
ccccc set grid
      call grid(xmin,xmax,ymin,ymax,dx,dy)
ccccc set initial data
      call init(q,qa)
ccccc set terminal time tf
      pi=4.d0*datan(1.d0)
      tf=2.d0*pi/pspd*dble(ncir)
      t =0.d0
      istp  = 0
      if (fstp .gt. 0.) then
       dprnt = tf/fstp
      else
       dprnt = 0.d0
      endif
      nprnt  = 0
      dnext  = dprnt
ccccc Ouput the initial data or rest the initial data.
      if (fstr .eq. 0) then
       call ioio(nprnt,t,q,0)
      else
       nprnt = fstr
       call ioio(nprnt,t,q,1)
       dnext = t+dprnt
      endif
ccccc Starting the time evolution
      do 1000 while ( t .lt. tf .and. istp .le. nstp)
ccccc compute the max velocity and get dt
      amax = 0.d0
      bmax = 0.d0
      do i = ibeg, iend
      do j = jbeg, jend
        ux   = q(i,j,2)/q(i,j,1)
        uy   = q(i,j,3)/q(i,j,1)
        amax = dmax1(amax,dabs(ux)+sspd)
        bmax = dmax1(bmax,dabs(uy)+sspd)
        if (amax<=(dabs(ux)+sspd))then
                dts%i = i
                dts%j = j
        endif
        if (bmax<=(dabs(ux)+sspd))then
                dts%i = i
                dts%j = j
        endif
      enddo
      enddo
      amax    = 1.1d0*amax
      bmax    = 1.1d0*bmax
      dt      = dmin1(dx/amax,dy/bmax)*dcfl
      if(dx/amax<=dy/bmax)then
        dts%direction = 0
      else
        dts%direction = 1
      endif
      dts%dt = dt
      dts%id = myid
ccccc MPI parallel computation require global mininum dt
      call LLSC(dt)
      call exchange_dts()
ccccc redefine dt if it is required
      if ( t + dt .gt. tf ) then
           dt = tf - t
           dts%dt = dt
      endif
      hdt=0.5d0*dt
ccccc MPI parallel computation change the buffer values
      call LLBUFFER(q)
ccccc compute boundary flux
      call bdry(t,q)
ccccc update flux1
      call flux1(dt,dx,dy,q,qa)
ccccc update source term
      call srce(t,dt,dx,dy,q,qa)
      t=t+dt
ccccc MPI parallel computation change the buffer values
      call LLBUFFER(qa)
ccccc compute boundary flux
      call bdry(t,qa)
ccccc update flux2
      call flux2(hdt,dx,dy,qa,q)
ccccc update source term
      call srce(t,hdt,dx,dy,qa,q)
ccccc Output result
       if ( fstp .eq. 0) then
        nprnt = 1 + nprnt
        call ioio(nprnt,t,q,0)
       elseif ( t.gt.dnext .or. t.ge.tf ) then
        nprnt = 1 + nprnt
        dnext = dnext + dprnt
        call ioio(nprnt,t,q,0)
       endif
       istp = istp + 1
ccccc On screen
      if (myid.eq.0) then
        write(6,*) 'istp = ', istp, ' t=',t, ' dt=', dt
        if(dts%dt<10.d-7)then
         write(*,'(a,g14.8,g14.8,g14.6,g1.1,g14.6)') ' dt find at (x,y)'
     c   ,x(dts%i),y(dts%j),'direct ',dts%direction,dts%dt
         write(*,*)

        endif
      endif
1000  continue
      CALL end_kepler
      stop
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  get variables get variables get variables get variables    
C  get variables get variables get variables get variables    
C  get variables get variables get variables get variables    
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine vars(xmin,xmax,ymin,ymax,
     &                dcfl,ncir,psrl,ichk,
     &                fstr,fstp,nstp,unk1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameters'
ccccc read variables from the archive variables
      open(3,file='variables')
      read(3,*)xmin,xmax,ymin,ymax,
     &         pspd,sspd,v0sp,bfam,
     &         dcfl,ncir,psrl,ichk,
     &         fstr,fstp,nstp,unk1 
      close(3)
ccccc reset bfam
ccccc set ph0
ccccc Huntley bar force phi and dphi
c        r    = 3.
cc       Op    = 120.
c       Op    =  pspd
c        a1   = 2.0d0
c       rsq   = r**2
cc       ar    = dsqrt(a1**2+rsq)
cc       aarr  = a1+ar
cc       p1    = 0.5d0*rsq/ar**3
cc       p2    = ( a1+2.d0*ar )/aarr**2
cc       psi   = p1*p2
c       psi   = rsq/(a1**2+rsq)**2
cc       dp1   = (r-(1.5d0*r**3)/ar**2)/ar**3
cc       dp2   = 2.d0*r*(1.d0-(a1+2.*ar)/aarr)/ar/(a1+ar)**2
cc       dpsi  = ((dp1*p2)+(p1*dp2))
c       dpsi  = 2.d0*r*(a1**2-r**2)/(a1**2+rsq)**3
ccccc Elmegreen rotation curve
c       rs    =  0.153d0
c       A~    = -0.081d0
c       B~    = -0.103d0
c       v0    =  230.4d0
c       rdrs  = r/rs
c       ycc   = v0/(rdrs**B+rdrs**(1.d0-A))/rs
c       Or    = ycc
c       ycc   = -r*dpsi+2.d0*psi*Or/(Op-Or)
c      bfam  = -bfam/ycc*r*Or**2
ccccc report by myid =0
!      if (myid.eq.0) then
ccccc report variables
!      write(6,*)'xmin=',xmin,' xmax=',xmax
!      write(6,*)'ymin=',ymin,' ymax=',ymax
!      write(6,*)'pspd=',pspd,' sspd=',sspd
!      write(6,*)'v0sp=',v0sp,' bfam=',bfam
!      write(6,*)'dcfl=',dcfl,' ncir=',ncir
!      write(6,*)'psrl=',psrl,' ichk=',ichk
!      write(6,*)'fstr=',fstr,' fstp=',fstp
!      write(6,*)'nstp=',nstp,' unk1=',unk1
!      endif
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  setup grid  setup grid  setup grid  setup grid  setup grid 
C  setup grid  setup grid  setup grid  setup grid  setup grid 
C  setup grid  setup grid  setup grid  setup grid  setup grid 
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine grid(xmin,xmax,ymin,ymax,dx,dy)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameters'
ccccc Cartesian grid type
ccccc global domain in y-direction
      ymax=dble(ymax)
      ymin=dble(ymin)
ccccc generate uniform mesh y for parallel/non-parallel computation
      ny =  jend-jbeg +1
      dy = (ymax-ymin)/dble(ny)
      do j = jbeg-jbuf-1,jend+jbuf
       yl(j) = ymin+dble(j)*dy
      enddo
      do j = jbeg-jbuf,jend+jbuf
        y(j) = 0.5d0*(yl(j-1)+yl(j))
      enddo
ccccc global domain in x-direction
      xmax=dble(xmax)
      xmin=dble(xmin)
ccccc MPI parallel computation:grid generation
ccccc uniform mesh in x-direction
      nx =  ieLL-ibeg +1
      dx = (xmax-xmin)/dble(nx)
      do i = ibeg-ibuf-1,ieLL+ibuf
        xxxl(i) = xmin+dble(i)*dx
      enddo
      do i = ibeg-ibuf,ieLL+ibuf
        xxx(i) = 0.5d0*(xxxl(i-1)+xxxl(i))
      enddo
      idx = (iend-ibeg+1)*myid
      do i= ibeg-ibuf-1,iend+ibuf
       xl(i) = xxxl(idx+i)
      enddo
      do i= ibeg-ibuf,iend+ibuf
        x(i) = xxx(idx+i)
      enddo
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  setup init  setup init  setup init  setup init  setup init 
C  setup init  setup init  setup init  setup init  setup init 
C  setup init  setup init  setup init  setup init  setup init 
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine init(q,qa)
      USE find_dt
      USE KEPLER,only: init_kepler
      USE rotation, only: angular
      USE density
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameters'
      dimension  q(ibeg-ibuf  :iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
      dimension qa(ibeg-ibuf  :iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
ccccc rcType(rc_def) is called to set the rotation curve
ccccc Elmegreen Rotation Curve for Cartesian Coordinates
ccccc variables in Elmegreen rotation curve rs,A,B,v0sp
ccccc get variables
      call vars(xmin,xmax,ymin,ymax,
     &          dcfl,ncir,psrl,ichk,
     &          fstr,fstp,nstp,unk1)

ccccc set nearly flat rotation curve
      call init_density('bell',sspd,10.d0)
      do i=ibeg-ibuf,iend+ibuf
      do j=jbeg-jbuf,jend+jbuf
             r   = sqrt(x(i)**2+y(j)**2)
             ycc = angular(r)
             rho = density_distribution(r)
              ux = -y(j)*ycc
              uy =  x(i)*ycc
       orsq(i,j) = ycc**2 + alpcsq2
       q(i,j,1)  = rho
       q(i,j,2)  = rho*ux
       q(i,j,3)  = rho*uy
      enddo
      enddo
ccccc copy q to qa
      do i=ibeg-ibuf,iend+ibuf
      do j=jbeg-jbuf,jend+jbuf
       qa(i,j,1) = q(i,j,1)
       qa(i,j,2) = q(i,j,2)
       qa(i,j,3) = q(i,j,3)
      enddo
      enddo

ccccc Init kepler mod
      CALL init_kepler(30.d0,-2.d0)
ccccc Init find_dt mod
      CALL find_dt_init()

      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C update flux1 update flux1 update flux1 update flux1 update  
C update flux1 update flux1 update flux1 update flux1 update  
C update flux1 update flux1 update flux1 update flux1 update  
C
C Godunov Godunov Godunov Godunov Godunov Godunov Godunov     
C Godunov Godunov Godunov Godunov Godunov Godunov Godunov     
C Godunov Godunov Godunov Godunov Godunov Godunov Godunov     
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine flux1(dt,dx,dy,q,qa)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameters'
      dimension    q(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
      dimension    qa(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)

      dimension    VX(ibeg-ibuf:iend+ibuf,1:3)
      dimension    VY(jbeg-jbuf:jend+jbuf,1:3)

      dimension   dV(1:3)
      dimension  VXL(ibeg-ibuf:iend+ibuf,1:3)
      dimension  VXR(ibeg-ibuf:iend+ibuf,1:3)
      dimension  VXM(ibeg-ibuf:iend+ibuf,1:3)
      dimension fhat(ibeg-ibuf:iend+ibuf,1:3)

      dimension  VYL(jbeg-jbuf:jend+jbuf,1:3)
      dimension  VYR(jbeg-jbuf:jend+jbuf,1:3)
      dimension  VYM(jbeg-jbuf:jend+jbuf,1:3)
      dimension ghat(jbeg-jbuf:jend+jbuf,1:3)
      fslop(aa,bb,cc)=dmax1(dsign(1.d0,(bb-aa)*(cc-bb)),0.d0)*
     &  dsign(1.d0,(cc-aa))*dmin1(dabs(0.5d0*(cc-aa)),
     &  dmin1(2.d0*dabs((bb-aa)),2.d0*dabs((cc-bb))))

      csq   = sspd**2
      dtddx = dt/dx
      dtddy = dt/dy

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC x-direction flux1 sweeping--start x loop
      do j = jbeg, jend
         do i = ibeg-ibuf,iend+ibuf
          VX(i,1) = q(i,j,1)
          VX(i,2) = q(i,j,2)/q(i,j,1)
          VX(i,3) = q(i,j,3)/q(i,j,1)
         enddo
         do i = ibeg-ibuf+1, iend+ibuf-1
          do k = 1, 3
           dV(k) = fslop(VX(i-1,k),VX(i,k),VX(i+1,k))
          enddo
          VXL(i   ,1) = VX(i,1) + 0.5d0*dV(1)
          VXL(i   ,2) = VX(i,2) + 0.5d0*dV(2)
          VXL(i   ,3) = VX(i,3) + 0.5d0*dV(3)
          VXR(i-1 ,1) = VX(i,1) - 0.5d0*dV(1)
          VXR(i-1 ,2) = VX(i,2) - 0.5d0*dV(2)
          VXR(i-1 ,3) = VX(i,3) - 0.5d0*dV(3)
         enddo
         if (myid.eq.0) then
         nrmibeg = ibeg
         do k = 1, 3
          dV(       k)= fslop(VX(ibeg,k),VX(ibeg+1,k),VX(ibeg+2,k))
          fhat(nrmibeg-1,k) = flxl(j,k)
         enddo
         VXL(ibeg  ,1)= VX(ibeg,1) + 0.5d0*dV(1)
         VXL(ibeg  ,2)= VX(ibeg,2) + 0.5d0*dV(2)
         VXL(ibeg  ,3)= VX(ibeg,3) + 0.5d0*dV(3)
         else
         nrmibeg = ibeg-1
         endif
         if (myid.eq.(nprocs-1)) then
         nrmiend = iend-1
         do k = 1, 3
          dV(       k) = fslop(VX(iend-2,k),VX(iend-1,k),VX(iend,k))
          fhat(nrmiend+1,k) = flxr(j,k)
         enddo
         VXR(iend-1,1)= VX(iend,1)-0.5d0*dV(1)
         VXR(iend-1,2)= VX(iend,2)-0.5d0*dV(2)
         VXR(iend-1,3)= VX(iend,3)-0.5d0*dV(3)
         else
         nrmiend = iend
         endif
         call xriem(ibeg,iend,ibuf,nrmibeg,nrmiend,sspd,VXL,VXR,VXM,j)
         do i = nrmibeg, nrmiend
          fhat(i  ,1) = VXM(i,1)*(          VXM(i,2)       )
          fhat(i  ,2) = VXM(i,1)*( VXM(i,2)*VXM(i,2) + csq )
          fhat(i  ,3) = VXM(i,1)*( VXM(i,3)*VXM(i,2)       )
         enddo
         do k = 1,3
          do i = ibeg, iend
           qa(i,j,k)=q(i,j,k)-dtddx*(fhat(i,k)-fhat(i-1,k))
          enddo
         enddo
      enddo
CCCCC x-direction flux1 sweeping--closing x loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC y-direction flux1 sweeping--starting y loop
      do i = ibeg, iend
         do j = jbeg-jbuf,jend+jbuf
            VY(j,1) = q(i,j,1)
            VY(j,2) = q(i,j,2)/q(i,j,1)
            VY(j,3) = q(i,j,3)/q(i,j,1)
         enddo
         do j = jbeg-jbuf+1, jend+jbuf-1
          do k = 1, 3
           dV(     k) = fslop(VY(j-1,k),VY(j,k),VY(j+1,k))
          enddo
          VYL(j   ,1) = VY(j,1) + 0.5d0*dV(1)
          VYL(j   ,2) = VY(j,2) + 0.5d0*dV(2)
          VYL(j   ,3) = VY(j,3) + 0.5d0*dV(3)
          VYR(j-1 ,1) = VY(j,1) - 0.5d0*dV(1)
          VYR(j-1 ,2) = VY(j,2) - 0.5d0*dV(2)
          VYR(j-1 ,3) = VY(j,3) - 0.5d0*dV(3)
         enddo
         nrmjbeg = jbeg
         do k = 1, 3
           dV(       k) = fslop(VY(jbeg,k),VY(jbeg+1,k),VY(jbeg+2,k))
           ghat(nrmjbeg-1,   k) = flxd(i,k);
         enddo
         VYL(jbeg  ,1) = VY(jbeg,1) + 0.5d0*dV(1)
         VYL(jbeg  ,2) = VY(jbeg,2) + 0.5d0*dV(2)
         VYL(jbeg  ,3) = VY(jbeg,3) + 0.5d0*dV(3)
         nrmjend = jend-1
         do k = 1, 3
           dV(       k) = fslop(VY(jend-2,k),VY(jend-1,k),VY(jend,k))
           ghat(nrmjend+1, k) = flxu(i,k);
         enddo
         VYR(jend-1,1) = VY(jend,1)-0.5d0*dV(1)
         VYR(jend-1,2) = VY(jend,2)-0.5d0*dV(2)
         VYR(jend-1,3) = VY(jend,3)-0.5d0*dV(3)
         call yriem(jbeg,jend,jbuf,nrmjbeg,nrmjend,sspd,VYL,VYR,VYM,i)
         do j = nrmjbeg, nrmjend
          ghat(j  ,1) = VYM(j,1)*(          VYM(j,3)       )
          ghat(j  ,2) = VYM(j,1)*( VYM(j,2)*VYM(j,3)       )
          ghat(j  ,3) = VYM(j,1)*( VYM(j,3)*VYM(j,3) + csq )
         enddo
         do k = 1,3
         do j = jbeg, jend
          qa(i,j,k) = qa(i,j,k) - dtddy*(ghat(j,k)-ghat(j-1,k))
         enddo
         enddo
      enddo
CCCCC y-direction flux1 sweeping--closing y loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C update flux2 update flux2 update flux2 update flux2 update  
C update flux2 update flux2 update flux2 update flux2 update  
C update flux2 update flux2 update flux2 update flux2 update  
C
C Godunov Godunov Godunov Godunov Godunov Godunov Godunov     
C Godunov Godunov Godunov Godunov Godunov Godunov Godunov     
C Godunov Godunov Godunov Godunov Godunov Godunov Godunov     
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine flux2(dt,dx,dy,q,qa)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameters'
      dimension    q(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
      dimension    qa(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)

      dimension    VX(ibeg-ibuf:iend+ibuf,1:3)
      dimension    VY(jbeg-jbuf:jend+jbuf,1:3)

      dimension   dV(1:3)
      dimension  VXL(ibeg-ibuf:iend+ibuf,1:3)
      dimension  VXR(ibeg-ibuf:iend+ibuf,1:3)
      dimension  VXM(ibeg-ibuf:iend+ibuf,1:3)
      dimension fhat(ibeg-ibuf:iend+ibuf,1:3)

      dimension  VYL(jbeg-jbuf:jend+jbuf,1:3)
      dimension  VYR(jbeg-jbuf:jend+jbuf,1:3)
      dimension  VYM(jbeg-jbuf:jend+jbuf,1:3)
      dimension ghat(jbeg-jbuf:jend+jbuf,1:3)
      fslop(aa,bb,cc)=dmax1(dsign(1.d0,(bb-aa)*(cc-bb)),0.d0)*
     &  dsign(1.d0,(cc-aa))*dmin1(dabs(0.5d0*(cc-aa)),
     &  dmin1(2.d0*dabs((bb-aa)),2.d0*dabs((cc-bb))))

      csq   = sspd**2
      dtddx = dt/dx
      dtddy = dt/dy

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC y-direction flux2 sweeping--starting y loop
      do i = ibeg, iend
         do j = jbeg-jbuf,jend+jbuf
            VY(j,1) = q(i,j,1)
            VY(j,2) = q(i,j,2)/q(i,j,1)
            VY(j,3) = q(i,j,3)/q(i,j,1)
         enddo
         do j = jbeg-jbuf+1, jend+jbuf-1
          do k = 1, 3
           dV(     k) = fslop(VY(j-1,k),VY(j,k),VY(j+1,k))
          enddo
          VYL(j   ,1) = VY(j,1) + 0.5d0*dV(1)
          VYL(j   ,2) = VY(j,2) + 0.5d0*dV(2)
          VYL(j   ,3) = VY(j,3) + 0.5d0*dV(3)
          VYR(j-1 ,1) = VY(j,1) - 0.5d0*dV(1)
          VYR(j-1 ,2) = VY(j,2) - 0.5d0*dV(2)
          VYR(j-1 ,3) = VY(j,3) - 0.5d0*dV(3)
         enddo
         nrmjbeg = jbeg
         do k = 1, 3
           dV(       k) = fslop(VY(jbeg,k),VY(jbeg+1,k),VY(jbeg+2,k))
           ghat(nrmjbeg-1,   k) = flxd(i,k);
         enddo
         VYL(jbeg  ,1) = VY(jbeg,1) + 0.5d0*dV(1)
         VYL(jbeg  ,2) = VY(jbeg,2) + 0.5d0*dV(2)
         VYL(jbeg  ,3) = VY(jbeg,3) + 0.5d0*dV(3)
         nrmjend = jend-1
         do k = 1, 3
           dV(       k) = fslop(VY(jend-2,k),VY(jend-1,k),VY(jend,k))
           ghat(nrmjend+1, k) = flxu(i,k);
         enddo
         VYR(jend-1,1) = VY(jend,1)-0.5d0*dV(1)
         VYR(jend-1,2) = VY(jend,2)-0.5d0*dV(2)
         VYR(jend-1,3) = VY(jend,3)-0.5d0*dV(3)
         call yriem(jbeg,jend,jbuf,nrmjbeg,nrmjend,sspd,VYL,VYR,VYM,i)
        do j = nrmjbeg, nrmjend
         ghat(j  ,1) = VYM(j,1)*(          VYM(j,3)       )
         ghat(j  ,2) = VYM(j,1)*( VYM(j,2)*VYM(j,3)       )
         ghat(j  ,3) = VYM(j,1)*( VYM(j,3)*VYM(j,3) + csq )
       enddo
       do k = 1,3
       do j = jbeg, jend
          qa(i,j,k) = 0.5*(qa(i,j,k)+q(i,j,k))
     &              - dtddy*(ghat(j,k)-ghat(j-1,k))
       enddo
       enddo
      enddo
CCCCC y-direction flux2 sweeping--closing y loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC x-direction flux2 sweeping--start x loop
      do j = jbeg, jend
         do i = ibeg-ibuf,iend+ibuf
          VX(i,1) = q(i,j,1)
          VX(i,2) = q(i,j,2)/q(i,j,1)
          VX(i,3) = q(i,j,3)/q(i,j,1)
         enddo
         do i = ibeg-ibuf+1, iend+ibuf-1
          do k = 1, 3
           dV(k) = fslop(VX(i-1,k),VX(i,k),VX(i+1,k))
          enddo
          VXL(i   ,1) = VX(i,1) + 0.5d0*dV(1)
          VXL(i   ,2) = VX(i,2) + 0.5d0*dV(2)
          VXL(i   ,3) = VX(i,3) + 0.5d0*dV(3)
          VXR(i-1 ,1) = VX(i,1) - 0.5d0*dV(1)
          VXR(i-1 ,2) = VX(i,2) - 0.5d0*dV(2)
          VXR(i-1 ,3) = VX(i,3) - 0.5d0*dV(3)
         enddo
         if (myid.eq.0) then
         nrmibeg = ibeg
         do k = 1, 3
          dV(       k)= fslop(VX(ibeg,k),VX(ibeg+1,k),VX(ibeg+2,k))
          fhat(nrmibeg-1,k) = flxl(j,k)
         enddo
         VXL(ibeg  ,1)= VX(ibeg,1) + 0.5d0*dV(1)
         VXL(ibeg  ,2)= VX(ibeg,2) + 0.5d0*dV(2)
         VXL(ibeg  ,3)= VX(ibeg,3) + 0.5d0*dV(3)
         else
         nrmibeg = ibeg-1
         endif
         if (myid.eq.(nprocs-1)) then
         nrmiend = iend-1
         do k = 1, 3
          dV(       k) = fslop(VX(iend-2,k),VX(iend-1,k),VX(iend,k))
          fhat(nrmiend+1,k) = flxr(j,k)
         enddo
         VXR(iend-1,1)= VX(iend,1)-0.5d0*dV(1)
         VXR(iend-1,2)= VX(iend,2)-0.5d0*dV(2)
         VXR(iend-1,3)= VX(iend,3)-0.5d0*dV(3)
         else
         nrmiend = iend
         endif
         call xriem(ibeg,iend,ibuf,nrmibeg,nrmiend,sspd,VXL,VXR,VXM,j)
         do i = nrmibeg, nrmiend
          fhat(i  ,1) = VXM(i,1)*(          VXM(i,2)       )
          fhat(i  ,2) = VXM(i,1)*( VXM(i,2)*VXM(i,2) + csq )
          fhat(i  ,3) = VXM(i,1)*( VXM(i,3)*VXM(i,2)       )
         enddo
         do k = 1,3
          do i = ibeg, iend
           qa(i,j,k)=qa(i,j,k)-dtddx*(fhat(i,k)-fhat(i-1,k))
          enddo
         enddo
      enddo
CCCCC x-direction flux2 sweeping--closing x loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  xriem solver xriem solver xriem solver xriem solver xriem  
C  xriem solver xriem solver xriem solver xriem solver xriem  
C  xriem solver xriem solver xriem solver xriem solver xriem  
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine xriem(ibeg,iend,ibuf,nbeg,nend,c,VL,VR,VM,j)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      integer   ibeg,iend,ibuf,nbeg,nend
      dimension VL(ibeg-ibuf:iend+ibuf,1:3)
      dimension VR(ibeg-ibuf:iend+ibuf,1:3)
      dimension VM(ibeg-ibuf:iend+ibuf,1:3)
      integer   n
      INTENT(IN)                                ::j
      INtr = 10

      do 1000 i = nbeg,nend
      if (VL(i,1).eq.VR(i,1).and.VL(i,2).eq.VR(i,2)) then
           VM(i,1) = VR(i,1)
           VM(i,2) = VR(i,2)
       if(VL(i,2).ge. 0.   ) then
           VM(i,3) = VL(i,3)
       else
           VM(i,3) = VR(i,3)
       endif
       goto 1000
      endif
      rhoL = VL(i,1)
      rhoR = VR(i,1)
       uxL = VL(i,2)
       uxR = VR(i,2)
       uyL = VL(i,3)
       uyR = VR(i,3)

      shkR = -c*(rhoR-rhoL)/dsqrt(rhoR*rhoL)
      rarR = -c*dlog(rhoR/rhoL)
       csq =  c*c

      if( rhoL .ge. rhoR ) then
        if ((uxR-uxL) .gt. rarR) then
          goto 100
        elseif ((uxR-uxL) .gt. -shkR) then
          goto 200
        else
          goto 400
        endif
      else
        if ((uxR-uxL) .lt. shkR ) then
          goto 400
        elseif ((uxR-uxL) .lt. -rarR) then
          goto 300
        else
          goto 100
        endif
      endif

cccccc 1-rarefaction wave, 2-rarefaction wave
cccccc 1-rarefaction wave, 2-rarefaction wave
cccccc 1-rarefaction wave, 2-rarefaction wave
100   ustar   = 0.5d0*(uxL+uxR)+0.5d0*c*dlog(rhoL/rhoR)
      rhostar = rhoL*dexp(-(ustar-uxL)/c)

      if ((uxL-c) .ge. 0.) then
         VM(i,2) = uxL
         VM(i,1) = rhoL
      elseif ((ustar-c) .ge. 0.) then
         VM(i,2) =  c
         VM(i,1) =  rhoL*dexp((uxL-c)/c)
      elseif ((ustar+c) .ge. 0.) then
         VM(i,2) =  ustar
         VM(i,1) =  rhostar
      elseif ((  uxR+c) .ge. 0.) then
         VM(i,2) =  -c
         VM(i,1) =  rhoR*dexp(-(uxR+c)/c)
      else
         VM(i,2) =   uxR
         VM(i,1) =  rhoR
      endif
      if ( VM(i,2) .ge. 0.) then
         VM(i,3) =  VL(i,3)
      else
         VM(i,3) =  VR(i,3)
      endif
      goto 1000
cccccc1-rarefaction wave, 2-shock wave
cccccc1-rarefaction wave, 2-shock wave
cccccc1-rarefaction wave, 2-shock wave
200   rho0 = rhoR
      do n = 1, INtr
        rhosR    = rho0/rhoR
        srhosR   = dsqrt(rhosR)
        aisrhosR = 1.d0/srhosR
        rhostar  = rho0 -
     &  (uxR-uxL+c*(srhosR-aisrhosR)+c*dlog(rho0/rhoL))
     &  /(0.5d0*c/rho0*(srhosR+aisrhosR)+c/rho0)
        if (dabs(rho0-rhostar) .lt. 1.e-12) then
            goto 210
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1R2S'
      write(6,*) 'rhostar,rhoL=',rshostar,rhoL
      write(6,*) 'uxL    , uxR=',uxL,      uxR
      call newton_exception(i,j,'x')
      stop
210   ustar  = uxL - c*dlog(rhostar/rhoL)
      sigma2 = uxR + c*dsqrt(rhostar/rhoR)
      if ( sigma2 .le. 0.) then
         VM(i,2) =  uxR
         VM(i,1) = rhoR
      elseif ((ustar-c) .le. 0.) then
         VM(i,2) = ustar
         VM(i,1) = rhostar
      elseif ((uxL-c) .le. 0.) then
         VM(i,2) =  c
         VM(i,1) =  rhoL*dexp((uxL-c)/c)
      else
         VM(i,2) = uxL
         VM(i,1) = rhoL
      endif
      if ( VM(i,2) .ge. 0.) then
         VM(i,3) = VL(i,3)
      else
         VM(i,3) = VR(i,3)
      endif
      goto 1000
cccccc 1-shock wave, 2-rarefaction wave
cccccc 1-shock wave, 2-rarefaction wave
cccccc 1-shock wave, 2-rarefaction wave
300   rho0 = rhoL
      do n = 1, INtr
        rhosL    = rho0/rhoL
        srhosL   = dsqrt(rhosL)
        aisrhosL = 1.d0/srhosL
        rhostar  = rho0-
     & (uxR-uxL+c*(srhosL-aisrhosL)+c*dlog(rho0/rhoR))
     &  /(0.5d0*c/rho0*(srhosL+aisrhosL)+c/rho0)
        if (dabs(rho0-rhostar) .lt. 1.e-12) then
            goto 310
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1S2R'
      write(6,*) 'rhostar,rhoL=',rshostar,rhoL
      write(6,*) 'uxL    , uxR=',uxL,      uxR
      call newton_exception(i,j,'x')
      stop
310   ustar  = uxR + c*dlog(rhostar/rhoR)
      sigma1 = uxL - c*dsqrt(rhostar/rhoL)
      if (sigma1 .ge. 0.) then
         VM(i,2) = uxL
         VM(i,1) = rhoL
      elseif ((ustar+c) .ge. 0.) then
         VM(i,2) = ustar
         VM(i,1) = rhostar
      elseif ((uxR+c) .ge. 0.) then
         VM(i,2) =  -c
         VM(i,1) =  rhoR*dexp(-(uxR+c)/c)
      else
         VM(i,2) = uxR
         VM(i,1) = rhoR
      endif
      if ( VM(i,2) .ge. 0.) then
         VM(i,3) = VL(i,3)
      else
         VM(i,3) = VR(i,3)
      endif
      goto 1000
c       1-shock wave, 2-shock wave
c       1-shock wave, 2-shock wave
c       1-shock wave, 2-shock wave
400   aa       = 1.d0/dsqrt(rhoL)+1.d0/dsqrt(rhoR)
      bb       =    dsqrt(rhoL)+   dsqrt(rhoR)
      uu       = uxR-uxL
      srhostar = (-uu+dsqrt(uu**2+4.d0*csq*aa*bb))/(2.d0*c*aa)
      rhostar  = srhostar*srhostar
      ustar    = uxL - c*(dsqrt(rhostar/rhoL)-dsqrt(rhoL/rhostar))
      sigma1   = uxL - c*dsqrt(rhostar/rhoL)
      sigma2   = uxR + c*dsqrt(rhostar/rhoR)
      if     (sigma1 .ge. 0.) then
       VM(i,2) = uxL
       VM(i,1) = rhoL
      elseif (sigma2 .ge. 0.) then
       VM(i,2) = ustar
       VM(i,1) = rhostar
      else
       VM(i,2) = uxR
       VM(i,1) = rhoR
      endif
         if  (VM(i,2) .ge. 0.) then
       VM(i,3) = VL(i,3)
      else
       VM(i,3) = VR(i,3)
      endif
      goto 1000
1000  continue
      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  yriem solver yriem solver yriem solver yriem solver yriem  
C  yriem solver yriem solver yriem solver yriem solver yriem  
C  yriem solver yriem solver yriem solver yriem solver yriem  
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine yriem(jbeg,jend,jbuf,nbeg,nend,c,VL,VR,VM,j)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      integer jbeg,jend,jbuf,nbeg,nend
      dimension VL(jbeg-jbuf:jend+jbuf,3)
      dimension VR(jbeg-jbuf:jend+jbuf,3)
      dimension VM(jbeg-jbuf:jend+jbuf,3)
      INTENT(IN)                                ::j
      integer  n, INtr

      INtr = 10

      do 1000 i = nbeg,nend
        if (VL(i,1).eq.VR(i,1).and.VL(i,3).eq.VR(i,3)) then
           VM(i,1) = VR(i,1)
           VM(i,3) = VR(i,3)
        if (VL(i,3).ge. 0.) then
           VM(i,2) = VL(i,2)
        else
           VM(i,2) = VR(i,2)
        endif
          goto 1000
        endif
        rhoL = VL(i,1)
        rhoR = VR(i,1)
         uxL = VL(i,2)
         uxR = VR(i,2)
         uyL = VL(i,3)
         uyR = VR(i,3)

        shkR = -c*(rhoR-rhoL)/dsqrt(rhoR*rhoL)
        rarR = -c*dlog(rhoR/rhoL)
         csq = c*c

        if (rhoL .ge. rhoR) then
          if ((uyR-uyL) .gt. rarR) then
              goto 100
          elseif ((uyR-uyL) .gt. -shkR) then
              goto 200
          else
              goto 400
          endif
        else
           if ((uyR-uyL) .lt. shkR ) then
               goto 400
           elseif ((uyR-uyL) .lt. -rarR) then
               goto 300
           else
               goto 100
           endif
        endif

cccccccc 1-rarefaction wave, 2-rarefaction wave
cccccccc 1-rarefaction wave, 2-rarefaction wave
cccccccc 1-rarefaction wave, 2-rarefaction wave
100     ustar = 0.5d0*(uyL+uyR)+0.5d0*c*dlog(rhoL/rhoR)
        rhostar = rhoL*dexp(-(ustar-uyL)/c)
        if ((uyL-c) .ge. 0.) then
           VM(i,3) = uyL
           VM(i,1) = rhoL
        elseif ((ustar-c) .ge. 0.) then
           VM(i,3) =  c
           VM(i,1) =  rhoL*dexp((uyL-c)/c)
        elseif ((ustar+c) .ge. 0.) then
           VM(i,3) = ustar
           VM(i,1) = rhostar
        elseif ((uyR+c) .ge. 0.) then
           VM(i,3) =  -c
           VM(i,1) =  rhoR*dexp(-(uyR+c)/c)
        else
           VM(i,3) = uyR
           VM(i,1) = rhoR
        endif
        if ( VM(i,3) .ge. 0.) then
           VM(i,2) = VL(i,2)
        else
           VM(i,2) = VR(i,2)
        endif
        goto 1000
c       1-rarefaction wave, 2-shock wave
200     rho0 = rhoR
        do n = 1, INtr
          rhosR   = rho0/rhoR
          srhosR  = dsqrt(rhosR)
          aisrhosR= 1.d0/srhosR
          rhostar = rho0 -
     &    (uyR-uyL+c*(srhosR-aisrhosR)+c*dlog(rho0/rhoL))
     &    /(0.5d0*c/rho0*(srhosR+aisrhosR)+c/rho0)
          if (dabs(rho0-rhostar) .lt. 1.e-12) then
            goto 210
          endif
          rho0 = rhostar
        enddo
        write(6,*) 'Newton Interation is divergent at y-dir 1R2S'
        write(6,*) 'rhostar,rhoL=',rshostar,rhoL
        write(6,*) 'uxL    , uxR=',uxL,      uxR
        call newton_exception(i,j,'y')
        stop
210     ustar = uyL-c*dlog(rhostar/rhoL)
        sigma2= uyR+c*dsqrt(rhostar/rhoR)
        if (sigma2 .le. 0.) then
           VM(i,3) = uyR
           VM(i,1) = rhoR
        elseif ((ustar-c) .le. 0.) then
           VM(i,3) = ustar
           VM(i,1) = rhostar
        elseif ((uyL-c) .le. 0.) then
           VM(i,3) =  c
           VM(i,1) =  rhoL*dexp((uyL-c)/c)
        else
           VM(i,3) = uyL
           VM(i,1) = rhoL
        endif
        if ( VM(i,3) .ge. 0.) then
           VM(i,2) = VL(i,2)
        else
           VM(i,2) = VR(i,2)
        endif
        goto 1000
c       1-shock wave, 2-rarefaction wave
300     rho0 = rhoL
        do n = 1, INtr
          rhosL    = rho0/rhoL
          srhosL   = dsqrt(rhosL)
          aisrhosL = 1.d0/srhosL
          rhostar  = rho0-
     &    (uyR-uyL+c*(srhosL-aisrhosL)+c*dlog(rho0/rhoR))/
     &    (0.5d0*c/rho0*(srhosL+aisrhosL)+c/rho0)
          if (dabs(rho0-rhostar) .lt. 1.e-12) then
             goto 310
          endif
          rho0 = rhostar
        enddo
        write(6,*) 'Newton Interation is divergent at y-dir 1S2R'
        write(6,*) 'rhostar,rhoL=',rshostar,rhoL
        write(6,*) 'uxL    , uxR=',uxL,      uxR
        call newton_exception(i,j,'y')
        stop
310      ustar = uyR + c*dlog(rhostar/rhoR)
        sigma1 = uyL - c*dsqrt(rhostar/rhoL)
        if (sigma1 .ge. 0.) then
           VM(i,3) = uyL
           VM(i,1) = rhoL
        elseif ((ustar+c) .ge. 0.) then
           VM(i,3) = ustar
           VM(i,1) = rhostar
        elseif ((uyR+c) .ge. 0.) then
           VM(i,3) =  -c
           VM(i,1) =  rhoR*dexp(-(uyR+c)/c)
        else
           VM(i,3) = uyR
           VM(i,1) = rhoR
        endif
        if ( VM(i,3) .ge. 0.) then
           VM(i,2) = VL(i,2)
        else
           VM(i,2) = VR(i,2)
        endif
        goto 1000
c       1-shock wave, 2-shock wave
400     aa      = 1.d0/dsqrt(rhoL)+1.d0/dsqrt(rhoR)
        bb      =    dsqrt(rhoL)+  dsqrt(rhoR)
        uu      = uyR-uyL
        srhostar= (-uu+dsqrt(uu**2+4.d0*csq*aa*bb))/(2.d0*c*aa)
        rhostar = srhostar*srhostar
        ustar   = uyL - c*(dsqrt(rhostar/rhoL)-dsqrt(rhoL/rhostar))
        sigma1  = uyL - c* dsqrt(rhostar/rhoL)
        sigma2  = uyR + c* dsqrt(rhostar/rhoR)
        if ( sigma1 .ge. 0.) then
         VM(i,3)= uyL
         VM(i,1)= rhoL
        elseif ( sigma2 .ge. 0.) then
         VM(i,3)= ustar
         VM(i,1)= rhostar
        else
         VM(i,3)= uyR
         VM(i,1)= rhoR
        endif
        if ( VM(i,3) .ge. 0.) then
         VM(i,2)= VL(i,2)
        else
         VM(i,2)= VR(i,2)
        endif
        goto 1000
1000   continue
       return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  compute boundary flux compute buondary flux compute boundary
C  compute boundary flux compute buondary flux compute boundary
C  compute boundary flux compute buondary flux compute boundary
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine bdry(t,q)
      USE rotation,only: angular
      USE density,only: density_distribution
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameters'
      dimension   q(ibeg-ibuf  :iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
ccccc get variables
ccccc Boundry Condition for Elmegreen Rotation Curve
      c     = sspd
      csq   = sspd**2
      dtddx = dt/dx
      dtddy = dt/dy

ccccc Left boundary.
      do j = jbeg-jbuf,jend+jbuf
         i   = ibeg
         r   = dsqrt(x(i)**2+y(j)**2)
         rho = density_distribution(r)
          ux = -y(j)*angular(r)
          uy =  x(i)*angular(r)
          dW = q(i,j,1)-rho
          dU = q(i,j,2)-rho*ux
          dV = q(i,j,3)-rho*uy
         i   = ibeg-1
         r   = dsqrt(xl(i)**2+y(j)**2)
         rho = density_distribution(r)
          ux = -y(j)*angular(r)
          uy =  xl(i)*angular(r)
          z1 =  ( ( ux+c )*dW-dU )/2.d0/c
          z2 =    (-uy   )*dW+dV
          z3 =  (-( ux-c )*dW+dU )/2.d0/c
         if ( (ux+c).gt.0.d0) z3 = 0.d0
         if ( (ux  ).gt.0.d0) z2 = 0.d0
         if ( (ux-c).gt.0.d0) z1 = 0.d0
         Wb      = rho   +(        z1+       z3    )
         Ub      = rho*ux+( (ux-c)*z1+(ux+c)*z3    )
         Vb      = rho*uy+( (uy  )*z1+(uy  )*z3+z2 )
         flxl(j,1)= Ub
         flxl(j,2)= Wb*( (Ub/Wb)**2+csq )
         flxl(j,3)= Vb*(  Ub/Wb         )
      enddo
ccccc Right boundary.
      do j = jbeg-jbuf,jend+jbuf
         i   = iend
         r   = dsqrt(x(i)**2+y(j)**2)
         rho = density_distribution(r)
          ux = -y(j)*angular(r)
          uy =  x(i)*angular(r)
          dW = q(i,j,1)-rho
          dU = q(i,j,2)-rho*ux
          dV = q(i,j,3)-rho*uy
         i   = iend
         r   = dsqrt(xl(i)**2+y(j)**2)
         rho = density_distribution(r)
          ux = -y(j)*angular(r)
          uy =  xl(i)*angular(r)
          z1 =  ( ( ux+c )*dW-dU )/2.d0/c
          z2 =    (-uy   )*dW+dV
          z3 =  (-( ux-c )*dW+dU )/2.d0/c
         if ( (ux+c).lt.0.d0) z3 = 0.d0
         if ( (ux  ).lt.0.d0) z2 = 0.d0
         if ( (ux-c).lt.0.d0) z1 = 0.d0
         Wb      = rho   +(        z1+       z3    )
         Ub      = rho*ux+( (ux-c)*z1+(ux+c)*z3    )
         Vb      = rho*uy+( (uy  )*z1+(uy  )*z3+z2 )
         flxr(j,1)= Ub
         flxr(j,2)= Wb*( (Ub/Wb)**2+csq )
         flxr(j,3)= Vb*(  Ub/Wb         )
      enddo
ccccc Down boundary.
      do i = ibeg-ibuf,iend+ibuf
         j   = jbeg
         r   = dsqrt(x(i)**2+y(j)**2)
         rho = density_distribution(r)
          ux = -y(j)*angular(r)
          uy =  x(i)*angular(r)
          dW = q(i,j,1)-rho
          dU = q(i,j,2)-rho*ux
          dV = q(i,j,3)-rho*uy
         j   = jbeg-1
         r   = dsqrt(x(i)**2+yl(j)**2)
         rho = density_distribution(r)
          ux = -yl(j)*angular(r)
          uy =  x(i)*angular(r)
          z1 =  ( ( uy+c )*dW-dV )/2.d0/c
          z2 =    (-ux   )*dW+dU
          z3 =  (-( uy-c )*dW+dV )/2.d0/c
         if ( (uy+c).gt.0.d0) z3 = 0.d0
         if ( (uy  ).gt.0.d0) z2 = 0.d0
         if ( (uy-c).gt.0.d0) z1 = 0.d0
         Wb      = rho   +(        z1+       z3    )
         Ub      = rho*ux+( (ux  )*z1+(ux  )*z3+z2 )
         Vb      = rho*uy+( (uy-c)*z1+(uy+c)*z3    )
         flxd(i,1)= Vb
         flxd(i,2)= Ub*(  Vb/Wb         )
         flxd(i,3)= Wb*( (Vb/Wb)**2+csq )
      enddo
ccccc Up boundary.
      do i = ibeg-ibuf,iend+ibuf
         j   = jend
         r   = dsqrt(x(i)**2+y(j)**2)
         rho = density_distribution(r)
          ux = -y(j)*angular(r)
          uy =  x(i)*angular(r)
          dW = q(i,j,1)-rho
          dU = q(i,j,2)-rho*ux
          dV = q(i,j,3)-rho*uy
         j   = jend
         r   = dsqrt(x(i)**2+yl(j)**2)
         rho = density_distribution(r)
          ux = -yl(j)*angular(r)
          uy =  x(i)*angular(r)
          z1 =  ( ( uy+c )*dW-dV )/2.d0/c
          z2 =    (-ux   )*dW+dU
          z3 =  (-( uy-c )*dW+dV )/2.d0/c
         if ( (uy+c).lt.0.d0) z3 = 0.d0
         if ( (uy  ).lt.0.d0) z2 = 0.d0
         if ( (uy-c).lt.0.d0) z1 = 0.d0
         Wb      = rho   +(        z1+       z3    )
         Ub      = rho*ux+( (ux  )*z1+(ux  )*z3+z2 )
         Vb      = rho*uy+( (uy-c)*z1+(uy+c)*z3    )
         flxu(i,1)= Vb
         flxu(i,2)= Ub*(  Vb/Wb         )
         flxu(i,3)= Wb*( (Vb/Wb)**2+csq )
      enddo
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  update source update source update source update source 
C  update source update source update source update source 
C  update source update source update source update source 
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine srce(t,dt,dx,dy,q,qa)
      USE force_from_companion_compoents
      USE kepler
      USE rotation, pspdf=>pspd
      USE smooth
      USE companion_force
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameters'
      dimension    q(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
      dimension    qa(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)

      REAL*8                    ::force(2)
      REAL*8                    ::mass_ratio

ccccc force contributed from star potential
      do i = ibeg,iend
      do j = jbeg,jend
        rho      = q(i,j,1)
        qa(i,j,2)= qa(i,j,2) - dt*rho*x(i)*orsq(i,j)
        qa(i,j,3)= qa(i,j,3) - dt*rho*y(j)*orsq(i,j)
      enddo
      enddo
ccccc force contributed from bar potential


ccccc Huntly bar potential for m=2 on Cartesian grid
c     bfam2 = 18224.0869
c     p02   = bfam2*dmin1(t/tf*10.d0,1.d0)
c     do i = ibeg,iend
c     do j = jbeg,jend
c         a2 = 5.2
c      rho   = q(i,j,1)
c     rsq    = x(i)**2+y(j)**2
c       r    = dsqrt(rsq)
c       p1   = rsq/(a2**2+rsq)**2
c      psi   = p02*p1
c      dpsi  = p02*2.d0*r*(a2**2-rsq)/(a2**2+rsq)**3
c      cs2s  = (x(i)**2-y(j)**2)/rsq
c      sn2s  = 2.d0*y(j)*x(i)/rsq
c      cs2t  = cos(2.d0*pspd*t)
c      sn2t  = sin(2.d0*pspd*t)
c      cs2st = (cs2s*cs2t+sn2s*sn2t)
c      sn2st = (sn2s*cs2t-cs2s*sn2t)
c      px    = (dpsi*cs2st)*x(i)/r+2.d0*psi*sn2st*y(j)/rsq
c      py    = (dpsi*cs2st)*y(j)/r-2.d0*psi*sn2st*x(i)/rsq

c      qa(i,j,2)= qa(i,j,2) - dt*rho*px
c      qa(i,j,3)= qa(i,j,3) - dt*rho*py
c      enddo
c      enddo

ccccc companion component from fourier potential
      mass_ratio= 0.1d0
      Am    = 30.d0
      do i = ibeg,iend
      do j = jbeg,jend
              r = dsqrt(x(i)**2+y(j)**2)
              force    =
     &        force_from_companion(x(i),y(j),t,Am,mass_ratio,'012')

              rho      = q(i,j,1)
              qa(i,j,2)= qa(i,j,2) + dt*rho*force(1)
              qa(i,j,3)= qa(i,j,3) + dt*rho*force(2)
      enddo
      enddo


ccccc REAL  companion force
c     do i = ibeg,iend
c     do j = jbeg,jend

c      rho   = q(i,j,1)
c      force = companion_direct_force(x(i),y(j),mass_ratio,ieLL)
c    &        *exp_smooth(r,18.d0)
c    &        *linear_increase(t,tf)
c      qa(i,j,2)= qa(i,j,2) + dt*rho*force(1)
c      qa(i,j,3)= qa(i,j,3) + dt*rho*force(2)

c     enddo
c     enddo



      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  MPI parallel computation to exchange data                  
C  MPI parallel computation to exchange data                  
C  MPI parallel computation to exchange data                  
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine LLBUFFER(q)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
ccccc MPI parallel computation
      include 'mpif.h'
      include 'parameters'
      dimension q(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
      integer  istatus(MPI_STATUS_SIZE)
      integer    ireqr,ireqs,ierr
ccccc  i-direction
       ilen = ibuf*(jend-jbeg+2*jbuf+1)*3
       itag1= 3
       if (myid.ne.(nprocs-1)) then
       do I=1,ibuf
       do J=jbeg-jbuf,jend+jbuf
       do K=1,3
               s1ibuf(I,J,K) = q(iend-I+1,J,K)
       enddo
       enddo
       enddo
c.......in ascc
c       write (*,*) "test"
       call MPI_ISend(s1ibuf,ilen,MPI_DOUBLE_PRECISION,myid+1,itag1,
     &          MPI_Comm_world,ireqs,ierr)
c.......in asiaa
c       call MPI_ISend(s1ibuf,ilen,MPI_DOUBLE_PRECISION,myid+1,itag1,
c     &          MPI_Comm_world,ierr)
       call MPI_Wait(ireqs,istatus,ierr)
       endif
       if (myid.ne.0) then
       call MPI_IRecv(r1ibuf,ilen,MPI_DOUBLE_PRECISION,myid-1,itag1,
     &          MPI_Comm_world,ireqr,ierr)
       call MPI_Wait(ireqr,istatus,ierr)
       do I=1,ibuf
       do J=jbeg-jbuf,jend+jbuf
       do K=1,3
               q(ibeg-I,J,K) = r1ibuf(I,J,K)
       enddo
       enddo
       enddo
       endif
 
       itag2 = 4
       if (myid.ne.0) then
       do I=1,ibuf
       do J=jbeg-jbuf,jend+jbuf
       do K=1,3
               s2ibuf(I,J,K) = q(ibeg+I-1,J,K)
       enddo
       enddo
       enddo
c......in ascc
       call MPI_ISend(s2ibuf,ilen,MPI_DOUBLE_PRECISION,myid-1,itag2,
     &          MPI_Comm_world,ireqs,ierr)
c......in asiaa
c       call MPI_ISend(s2ibuf,ilen,MPI_DOUBLE_PRECISION,myid-1,itag2,
c     &          MPI_Comm_world,ierr)
       call MPI_Wait(ireqs,istatus,ierr)
       endif
       if (myid.ne.(nprocs-1)) then
       call MPI_IRecv(r2ibuf,ilen,MPI_DOUBLE_PRECISION,myid+1,itag2,
     &          MPI_Comm_world,ireqr,ierr)
       call MPI_Wait(ireqr,istatus,ierr)
       do I=1,ibuf
       do J=jbeg-jbuf,jend+jbuf
       do K=1,3
               q(iend+I,J,K) = r2ibuf(I,J,K)
       enddo
       enddo
       enddo
       endif
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  MPI parallel computation to find minimal dt among all nodes
C  MPI parallel computation to find minimal dt among all nodes
C  MPI parallel computation to find minimal dt among all nodes
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine LLSC(dt)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
ccccc MPI parallel computation
      include 'mpif.h'
      include 'parameters'
      dimension  sdt(0:0), rdt(0:0)
      dimension s0dt(0:0),r0dt(0:nprocs-1)
      integer  istatus(MPI_STATUS_SIZE)
      sdt(0)= dt
      itag1 = 1
      if (myid.ne.0) then
              call MPI_Send(sdt,1,MPI_DOUBLE_PRECISION,0,itag1,
     &        MPI_COMM_WORLD,ierr)
      endif
      if (myid.eq.0) then
              r0dt(0)=dt
              do I=1,nprocs-1
              call MPI_Recv(r0dt(I),1,MPI_DOUBLE_PRECISION,I,itag1,
     &        MPI_COMM_WORLD,istatus,ierr)
c              write(6,*) 'r0dt=',r0dt(I), ' myid=',I
              enddo
              dtmin = r0dt(0)
c              write(6,*) 'dtmin=',dtmin, ' myid=',myid
              do I=1,nprocs-1
               dtmin= dmin1(dtmin,r0dt(I))
              enddo
                dt    = dtmin
      endif

      !Update dt to everyone
      call MPI_Bcast(dt,1,MPI_DOUBLE_PRECISION,0,
     c               MPI_COMM_WORLD,ierr)
      !dt = rdt(0)

      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  ioio ioio ioio ioio ioio ioio ioio ioio ioio ioio ioio ioio
C  ioio ioio ioio ioio ioio ioio ioio ioio ioio ioio ioio ioio
C  ioio ioio ioio ioio ioio ioio ioio ioio ioio ioio ioio ioio
C
C2345678921234567893123456789412345678951234567896123456789712
      subroutine ioio(nprnt,t,q,IOC)
      USE kepler
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameters'
      include 'mpif.h'
      integer iiistatus(MPI_STATUS_SIZE)
      dimension q(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf,1:3)
      dimension tm(1)
      integer nprnt, IOC
      integer nx,ny,npt
      integer*8 matOpen, mxCreateFull,matClose
      integer*8 mxGetPr, matPutMatrix,matGetNextMatrix
      integer*8 qpt, tpt, spt, fpt
      integer ic1, ic2, ic3, ic4
      character*14 flnm
      tm(1)       = t
      nx          = (iend-ibeg+1+2*ibuf)
      ny          = (jend-jbeg+1+2*jbuf)
      npt         = 3*nx*ny
      ic1         = (nprnt          -mod(nprnt,1000))/1000
      ic2         = (mod(nprnt,1000)-mod(nprnt, 100))/100
      ic3         = (mod(nprnt, 100)-mod(nprnt,  10))/10
      ic4         =  mod(nprnt,  10)
      flnm( 1:6 ) = 'mats/M'
      flnm( 7:7 ) = char(ichar('0')+ic1)
      flnm( 8:8 ) = char(ichar('0')+ic2)
      flnm( 9:9 ) = char(ichar('0')+ic3)
      flnm(10:10) = char(ichar('0')+ic4)
      flnm(11:14)  = '.mat'
      itag1 = 5
      if (IOC .eq.0) then
      if (myid.ne.0) then
c......in ascc
c      call MPI_Send(q,3*nx*ny,MPI_DOUBLE_PRECISION,0,itag1,
c     &     MPI_COMM_WORLD,istatus,ierr)
c......in asiaa
      call MPI_Send(q,3*nx*ny,MPI_DOUBLE_PRECISION,0,itag1,
     &     MPI_COMM_WORLD,ierr)
      endif
      if (myid.eq.0) then
       do k = 1,3
       do i = ibeg-ibuf,iend+ibuf
       do j = jbeg-jbuf,jend+jbuf
           qLL(i,j,k) = q(i,j,k)
       enddo
       enddo
       enddo
       do n = 1,nprocs-1
           call MPI_Recv(qm,3*nx*ny,MPI_DOUBLE_PRECISION,n,itag1,
     &     MPI_COMM_WORLD,iiistatus,ierr)
           do k = 1,3
           do i = ibeg,iend+ibuf
           do j = jbeg-jbuf,jend+jbuf
           qLL((iend-ibeg+1)*n+i,j,k)=qm(i,j,k)
           enddo
           enddo
           enddo
       enddo
      nx          = (ieLL-ibeg+1+2*ibuf)
      ny          = (jend-jbeg+1+2*jbuf)
      npt         = 3*nx*ny
       fpt = matOpen(flnm,"w")
       qpt = mxCreateFull(1,npt,0)
       tpt = mxCreateFull(1,1,0)
       call  mxCopyReal8ToPtr(qLL,mxGetPr(qpt),npt)
       call  mxCopyReal8ToPtr(tm,mxGetPr(tpt),1)
       call  mxSetName(qpt,'q')
       call  mxSetName(tpt,'t')
       spt = matPutMatrix(fpt,qpt)
       spt = matPutMatrix(fpt,tpt)
        spt = matClose(fpt)
        call  mxFreeMatrix(qpt)
        call  mxFreeMatrix(tpt)
        call  print_kepler(t)
        write(*,*)'printing out frame:',nprnt
       endif
      else
       !read from file
       if (myid.eq.0) then
       nx  = (ieLL-ibeg+1+2*ibuf)
       ny  = (jend-jbeg+1+2*jbuf)
       npt = 3*nx*ny
       fpt = matOpen(flnm,'r')
       if (fpt .eq. 0) then
         write(6,*) 'To open ', flnm, ' fail....'
         stop
       end if
       qpt = matGetNextMatrix(fpt, 'q')
       tpt = matGetNextMatrix(fpt, 't' )
       call  mxCopyPtrToReal8(mxGetPr(qpt),qLL,npt)
       call  mxCopyPtrToReal8(mxGetPr(tpt),tm,1)
       t = tm(1)
       do k = 1,3
       do i = ibeg-ibuf,iend+ibuf
       do j = jbeg-jbuf,jend+jbuf
           q(i,j,k)=qLL(i,j,k)
       enddo
       enddo
       enddo
       nx          = (iend-ibeg+1+2*ibuf)
       ny          = (jend-jbeg+1+2*jbuf)
       npt         = 3*nx*ny
       do n = 1,nprocs-1
           do k = 1,3
           do i = ibeg,iend+ibuf
           do j = jbeg-jbuf,jend+jbuf
           qm(i,j,k)=qLL((iend-ibeg+1)*n+i,j,k)
           enddo
           enddo
           enddo
           call MPI_Send(qm,3*nx*ny,MPI_DOUBLE_PRECISION,n,itag1,
     &     MPI_COMM_WORLD,ierr)
           call MPI_Send(tm,1,MPI_DOUBLE_PRECISION,n,itag1,
     &     MPI_COMM_WORLD,ierr)
       enddo
       else
       nx          = (iend-ibeg+1+2*ibuf)
       ny          = (jend-jbeg+1+2*jbuf)
       npt         = 3*nx*ny
          call MPI_Recv(q,3*nx*ny,MPI_DOUBLE_PRECISION,0,itag1,
     &     MPI_COMM_WORLD,iiistatus,ierr)
          call MPI_Recv(tm,1,MPI_DOUBLE_PRECISION,0,itag1,
     &     MPI_COMM_WORLD,iiistatus,ierr)
       t = tm(1)

      endif
      endif
      end
