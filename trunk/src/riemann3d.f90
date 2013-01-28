subroutine riemann3d(q_loc,temp1_loc,temp2_loc,dd)
use common_params
implicit none
integer::dd
double precision::    q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR)
double precision::temp1_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR)
double precision::temp2_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR)
double precision::dtddx  ! dt/dx
double precision::fslop,vslop !slope limiter
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::vzL,vzM,vzR
double precision::eneL,eneM,eneR
double precision::pL,pM,pR
double precision,dimension(:,:,:,:),allocatable ::SL
double precision,dimension(:,:,:,:),allocatable ::Fx
double precision::Ftemp(NVAR), ql(NVAR), qr(NVAR), slopeL(NVAR), slopeR(NVAR)
integer::i,j,k
integer::cx,cy,cz,f1,f2,f3,f4,f5
double precision:: coef(NVAR),g1,g2,g3,g4,g5
#ifdef MHD
integer::f6,f7,f8,f9,f10,f11
double precision::g6,g7,g8,g9,g10,g11
double precision::bxL,bxM,bxR
double precision::byL,byM,byR
double precision::bzL,bzM,bzR
double precision,dimension(:,:,:),allocatable::Ex,Ey,Ez
#endif

allocate(SL(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR))
allocate(Fx(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR))

!!!!!!  dd == 1 ==> x, dd==2 ==> y, dd==3 ==> z  cccccc
      if (dd .eq. 1) then
         cx = 1
         cy = 0
         cz = 0
         f1 = 1
         f2 = 2
         f3 = 3
         f4 = 4
         f5 = 5
         coef(1) = 1.d0
         coef(2) = 1.d0
         coef(3) = 1.d0
         coef(4) = 1.d0
         coef(5) = 1.d0
         g1 = 1.d0
         g2 = 1.d0
         g3 = 1.d0
         g4 = 1.d0
         g5 = 1.d0
#ifdef MHD
         f5 = 5
         f6 = 6
         f7 = 7
         f8 = 8
         f9 = 9
        f10 = 10
        f11 = 11 
         coef(5) = 1.d0
         coef(6) = 1.d0
         coef(7) = 1.d0
         coef(8) = 1.d0
         g5 = 1.d0
         g6 = 1.d0
         g7 = 1.d0
         g8 = 1.d0
         g9 = 1.d0
        g10 = 1.d0
        g11 = 1.d0
#endif
      elseif (dd .eq. 2) then
         cx = 0
         cy = 1
         cz = 0
         f1 = 1
         f2 = 3
         f3 = 2
         f4 = 4
         f5 = 5
         coef(1) = 1.d0
         coef(2) = -1.d0
         coef(3) = 1.d0
         coef(4) = 1.d0
         coef(5) = 1.d0
         g1 = 1.d0
         g2 = 1.d0
         g3 = -1.d0
         g4 = 1.d0
         g5 = 1.d0
#ifdef MHD
         f5 = 6
         f6 = 5
         f7 = 7
         f8 = 8
         f9 = 10
        f10 = 9
        f11 = 11
         coef(5) = -1.d0
         coef(6) = 1.d0
         coef(7) = 1.d0
         coef(8) = 1.d0
         g5=1.d0
         g6=-1.d0
         g7=1.d0
         g8=1.d0
         g9=1.d0
        g10=-1.d0
        g11=1.d0
#endif
      elseif (dd .eq. 3) then
         cx=0
         cy=0
         cz=1
         f1=1
         f2=4
         f3=3
         f4=2
         f5=5
         coef(1)=1.d0
         coef(2)=-1.d0
         coef(3)=1.d0
         coef(4)=1.d0
         coef(5)=1.d0
         g1=1.d0
         g2=1.d0
         g3=1.d0
         g4=-1.d0
         g5=1.d0
#ifdef MHD
         f5=7
         f6=6
         f7=5
         f8=8
         f9=11
        f10=10
        f11=9
        coef(5)=-1.d0
        coef(6)=1.d0
        coef(7)=1.d0
        coef(8)=1.d0
         g5=1.d0
         g6=1.d0
         g7=-1.d0
         g8=1.d0
         g9=1.d0
        g10=1.d0
        g11=-1.d0
#endif
      endif
!cccccccccccc slope limiter  cccccccccccccccc
       do k=0, ncell_loc(3)+1
         do j=0, ncell_loc(2)+1
           do i=0, ncell_loc(1)+1
           rhoL = g1*q_loc(i-cx,j-cy,k-cz,f1)
           rhoM = g1*q_loc(i   ,j   ,k,   f1)
           rhoR = g1*q_loc(i+cx,j+cy,k+cz,f1)

           vxL = g2*q_loc(i-cx,j-cy,k-cz,f2)/rhoL
           vxM = g2*q_loc(i   ,j   ,k,   f2)/rhoM
           vxR = g2*q_loc(i+cx,j+cy,k+cz,f2)/rhoR

           vyL = g3*q_loc(i-cx,j-cy,k-cz,f3)/rhoL
           vyM = g3*q_loc(i   ,j   ,k,   f3)/rhoM
           vyR = g3*q_loc(i+cx,j+cy,k+cz,f3)/rhoR

           vzL = g4*q_loc(i-cx,j-cy,k-cz,f4)/rhoL
           vzM = g4*q_loc(i   ,j   ,k,   f4)/rhoM
           vzR = g4*q_loc(i+cx,j+cy,k+cz,f4)/rhoR

           eneL  = g5*q_loc(i-cx,j-cy,k-cz,f5)/rhoL
           eneM  = g5*q_loc(i   ,j   ,k,   f5)/rhoM
           eneR  = g5*q_loc(i+cx,j+cy,k+cz,f5)/rhoR

             pL = (rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2))*(gam-1.d0)
             pM = (rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2+vzM**2))*(gam-1.d0)
             pR = (rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2))*(gam-1.d0)
           
#ifdef MHD
            bxL = 0.5d0*(g5*q_loc(i-cx,j-cy,k-cz,f5)+g9*q_loc(i-cx,j-cy,k-cz,f9))
            bxM = 0.5d0*(g5*q_loc(i   ,j   ,k,   f5)+g9*q_loc(i   ,j   ,k   ,f9))
            bxR = 0.5d0*(g5*q_loc(i+cx,j+cy,k+cz,f5)+g9*q_loc(i+cx,j+cy,k+cz,f9))
 
            byL = 0.5d0*(g6*q_loc(i-cx,j-cy,k-cz,f6)+g10*q_loc(i-cx,j-cy,k-cz,f10))
            byM = 0.5d0*(g6*q_loc(i   ,j   ,k,   f6)+g10*q_loc(i   ,j   ,k   ,f10))
            byR = 0.5d0*(g6*q_loc(i+cx,j+cy,k+cz,f6)+g10*q_loc(i+cx,j+cy,k+cz,f10))

            bzL = 0.5d0*(g7*q_loc(i-cx,j-cy,k-cz,f7)+g11*q_loc(i-cx,j-cy,k-cz,f11))
            bzM = 0.5d0*(g7*q_loc(i   ,j   ,k,   f7)+g11*q_loc(i   ,j   ,k   ,f11))
            bzR = 0.5d0*(g7*q_loc(i+cx,j+cy,k+cz,f7)+g11*q_loc(i+cx,j+cy,k+cz,f11))

           eneL = g8*q_loc(i-cx,j-cy,k-cz,f8)/rhoL
           eneM = g8*q_loc(i   ,j   ,k,   f8)/rhoM
           eneR = g8*q_loc(i+cx,j+cy,k+cz,f8)/rhoR

             pL = (rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2)-0.5d0*(bxL**2+byL**2+bzL**2))*(gam-1.d0)
             pM = (rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2+vzM**2)-0.5d0*(bxM**2+byM**2+bzM**2))*(gam-1.d0)
             pR = (rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2)-0.5d0*(bxR**2+byR**2+bzR**2))*(gam-1.d0)
#endif


         SL(i,j,k,1) = vslop(rhoL,rhoM,rhoR)
         SL(i,j,k,2) = vslop(vxL,vxM,vxR)
         SL(i,j,k,3) = vslop(vyL,vyM,vyR)
         SL(i,j,k,4) = vslop(vzL,vzM,vzR)
         SL(i,j,k,5) = vslop(pL,pM,pR)
#ifdef MHD
         SL(i,j,k,5) = vslop( bxL, bxM, bxR)
         SL(i,j,k,6) = vslop( byL, byM, byR)
         SL(i,j,k,7) = vslop( bzL, bzM, bzR)
         SL(i,j,k,8) = vslop(pL,pM,pR)
#endif
           enddo
         enddo
       enddo
!cccccccccccc flux calculation cccccccccccccc
       do k=0, ncell_loc(3)+1
         do j=0, ncell_loc(2)+1
           do i=0, ncell_loc(1)+1

                   ql(1) = g1*q_loc(i,j,k,f1)
                   ql(2) = g2*q_loc(i,j,k,f2)/ql(1)
                   ql(3) = g3*q_loc(i,j,k,f3)/ql(1)
                   ql(4) = g4*q_loc(i,j,k,f4)/ql(1)
                   ql(5) = g5*(q_loc(i,j,k,f5)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2+ql(4)**2))*(gam-1.d0)

                   qr(1) = g1*q_loc(i+cx,j+cy,k+cz,f1)
                   qr(2) = g2*q_loc(i+cx,j+cy,k+cz,f2)/qr(1)
                   qr(3) = g3*q_loc(i+cx,j+cy,k+cz,f3)/qr(1)
                   qr(4) = g4*q_loc(i+cx,j+cy,k+cz,f4)/qr(1)
                   qr(5) = g5*(q_loc(i+cx,j+cy,k+cz,f5)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2+qr(4)**2))*(gam-1.d0)

#ifdef MHD
                   ql(5) = 0.5d0*(g5*q_loc(i,j,k,f5)+ g9*q_loc(i,j,k,f9 ))  ! bx
                   ql(6) = 0.5d0*(g6*q_loc(i,j,k,f6)+g10*q_loc(i,j,k,f10))  ! by
                   ql(7) = 0.5d0*(g7*q_loc(i,j,k,f7)+g11*q_loc(i,j,k,f11))  ! bz
                   ql(8) = g8*(q_loc(i,j,k,f8)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2+ql(4)**2)-0.5d0*(ql(5)**2+ql(6)**2+ql(7)**2))*(gam-1.d0)

                   qr(5) = 0.5d0*(g5*q_loc(i+cx,j+cy,k+cz,f5)+ g9*q_loc(i+cx,j+cy,k+cz,f9 ))
                   qr(6) = 0.5d0*(g6*q_loc(i+cx,j+cy,k+cz,f6)+g10*q_loc(i+cx,j+cy,k+cz,f10))
                   qr(7) = 0.5d0*(g7*q_loc(i+cx,j+cy,k+cz,f7)+g11*q_loc(i+cx,j+cy,k+cz,f11))
                   qr(8) = g8*(q_loc(i+cx,j+cy,k+cz,f8)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2+qr(4)**2)-0.5d0*(qr(5)**2+qr(6)**2+qr(7)**2))*(gam-1.d0)                   
#endif

               slopeL(:) = SL(i,j,k,:)
               slopeR(:) = SL(i+cx,j+cy,k+cz,:)
!!!!!! ccccccc puting if-endif into loop is very bad
!!!!!! try to use preprocessor instead !!!!
!           if(solver_type .eq. 1) then
!              call fluxRieIso2d(ql,qr,slopeL,slopeR,Ftemp)
!           elseif(solver_type .eq. 2) then
!              call fluxHLLC2d(ql,qr,slopeL,slopeR,Ftemp)
!           endif
#ifdef ISOTHERMAL
#ifdef ExactRie
  call fluxRieIso3d(ql,qr,slopeL,slopeR,Ftemp)
#endif
#ifdef HLL
  call fluxHLLIso3d(ql,qr,slopeL,slopeR,Ftemp)
#endif
#endif

#ifdef ADIABATIC
#ifdef HLLC
  call fluxHLLC3d(ql,qr,slopeL,slopeR,Ftemp)
#endif
#ifdef HLL
#ifdef MHD
  call fluxHLLAdiMHD(ql,qr,slopeL,slopeR,Ftemp)
#else
  call fluxHLLAdi3d(ql,qr,slopeL,slopeR,Ftemp)
#endif
#endif
#ifdef HLLD
  call fluxHLLDAdiMHD(ql,qr,slopeL,slopeR,Ftemp)
#endif
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 Fx(i,j,k,:) = Ftemp(:)
          enddo
        enddo
      enddo

       call bndflux3d(q_loc,Fx,dd)  ! implement your boundary flux in this subroutine

        dtddx = dt/dx
        do k=1, ncell_loc(3)
          do j=1, ncell_loc(2)
            do i=1, ncell_loc(1)
               temp2_loc(i,j,k,1) = temp1_loc(i,j,k,1)+ &
                 coef(1)*(Fx(i-cx,j-cy,k-cz,f1)-Fx(i,j,k,f1))*dtddx
               temp2_loc(i,j,k,2) = temp1_loc(i,j,k,2)+ &
                 coef(2)*(Fx(i-cx,j-cy,k-cz,f2)-Fx(i,j,k,f2))*dtddx
               temp2_loc(i,j,k,3) = temp1_loc(i,j,k,3)+ &
                 coef(3)*(Fx(i-cx,j-cy,k-cz,f3)-Fx(i,j,k,f3))*dtddx
               temp2_loc(i,j,k,4) = temp1_loc(i,j,k,4)+ &
                 coef(4)*(Fx(i-cx,j-cy,k-cz,f4)-Fx(i,j,k,f4))*dtddx
#ifdef MHD
               !temp2_loc(i,j,k,6) = temp1_loc(i,j,k,6)+ &
               !  coef(6)*(Fx(i-cx,j-cy,k-cz,f6)-Fx(i,j,k,f6))*dtddx
               !temp2_loc(i,j,k,7) = temp1_loc(i,j,k,7)+ &
               !  coef(7)*(Fx(i-cx,j-cy,k-cz,f7)-Fx(i,j,k,f7))*dtddx
               temp2_loc(i,j,k,8) = temp1_loc(i,j,k,8)+ &
                 coef(8)*(Fx(i-cx,j-cy,k-cz,f8)-Fx(i,j,k,f8))*dtddx
#else
               temp2_loc(i,j,k,5) = temp1_loc(i,j,k,5)+ &
                 coef(5)*(Fx(i-cx,j-cy,k-cz,f5)-Fx(i,j,k,f5))*dtddx
#endif
             enddo
          enddo
        enddo
#ifdef MHD
      allocate(Ex(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf))
      allocate(Ey(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf))
      allocate(Ez(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf))
        Ex=0.d0
        Ey=0.d0
        Ez=0.d0

      do k=0,ncell_loc(3)
        do j=0,ncell_loc(2)
          do i=0,ncell_loc(1)
             if(dd .eq. 1) then
               Ez(i,j,k)=0.25d0*(-Fx(i,j,k,f6)-Fx(i,j+1,k,f6))*coef(6)
               Ey(i,j,k)=0.25d0*( Fx(i,j,k,f7)+Fx(i,j,k+1,f7))*coef(7)
             endif

             if(dd .eq. 2) then
               Ez(i,j,k)=0.25d0*( Fx(i+1,j,k,f5)+Fx(i,j,k,f5))*coef(5)
               Ex(i,j,k)=0.25d0*(-Fx(i,j,k+1,f7)-Fx(i,j,k,f7))*coef(7)
             endif
  
             if(dd .eq. 3) then
               Ey(i,j,k)=0.25d0*(-Fx(i+1,j,k,f5)-Fx(i,j,k,f5))*coef(5)
               Ex(i,j,k)=0.25d0*( Fx(i,j+1,k,f6)+Fx(i,j,k,f6))*coef(6)
             endif
          enddo
        enddo
      enddo


     do k=1,ncell_loc(3)
       do j=1,ncell_loc(2)
         do i=1,ncell_loc(1)

       temp2_loc(i,j,k,5)  = temp1_loc(i,j,k,5)+ &
            (Ez(i-1,j-1,k)-Ez(i-1,j,k  ))*dtddx & 
            +(Ey(i-1,j  ,k)-Ey(i-1,j,k-1))*dtddx
       temp2_loc(i,j,k,9)  = temp1_loc(i,j,k,9)+ &
            (Ez(i  ,j-1,k)-Ez(i  ,j,k  ))*dtddx & 
            +(Ey(i  ,j  ,k)-Ey(i  ,j,k-1))*dtddx


       temp2_loc(i,j,k,6)=temp1_loc(i,j,k,6)+ &
            (Ez(i,j-1,k)-Ez(i-1,j-1,k))*dtddx &
            +(Ex(i,j-1,k-1)-Ex(i,j-1,k))*dtddx
       temp2_loc(i,j,k,10) = temp1_loc(i,j,k,10)+ &
            (Ez(i,j,k  )-Ez(i-1,j,k  ))*dtddx &
            +(Ex(i,j  ,k-1)-Ex(i,j  ,k))*dtddx

       temp2_loc(i,j,k,7)=temp1_loc(i,j,k,7)+ &
            (Ey(i-1,j,k-1)-Ey(i,j,k-1))*dtddx+ &
            (Ex(i,j,k-1)-Ex(i,j-1,k-1))*dtddx
       temp2_loc(i,j,k,11)=temp1_loc(i,j,k,11)+ &
            (Ey(i-1,j,k)-Ey(i,j,k))*dtddx+ &
            (Ex(i,j,k)-Ex(i,j-1,k))*dtddx


         enddo
       enddo
     enddo
    
      deallocate(Ex)
      deallocate(Ey)
      deallocate(Ez)
#endif

deallocate(SL)
deallocate(Fx)

end subroutine
