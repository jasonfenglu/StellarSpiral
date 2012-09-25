subroutine riemann2d(q_loc,temp1_loc,temp2_loc,dd)
use common_params
implicit none
integer::dd
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::temp1_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::temp2_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::dtddx  ! dt/dx
double precision::fslop,vslop !slope limiter
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::eneL,eneM,eneR
double precision::pL,pM,pR
double precision,dimension(:,:,:),allocatable ::SL
double precision,dimension(:,:,:),allocatable ::Fx
double precision::Ftemp(NVAR), ql(NVAR), qr(NVAR), slopeL(NVAR), slopeR(NVAR)
integer::i,j
integer::cx,cy,f1,f2,f3,f4
double precision:: coef(NVAR),g1,g2,g3,g4
#ifdef MHD
integer::f5,f6,f7,f8,f9,f10,f11
double precision::g5,g6,g7,g8,g9,g10,g11
double precision::vzL,vzM,vzR
double precision::bxL,bxM,bxR
double precision::byL,byM,byR
double precision::bzL,bzM,bzR
double precision,dimension(:,:),allocatable::Ex,Ey,Ez
#endif

allocate(SL(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR))
allocate(Fx(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR))

!!!!!!  dd == 1 ==> x, dd==2 ==> y  cccccc
      if (dd .eq. 1) then
         cx = 1
         cy = 0
         f1 = 1
         f2 = 2
         f3 = 3
         f4 = 4
         coef(1) = 1.d0
         coef(2) = 1.d0
         coef(3) = 1.d0
         coef(4) = 1.d0
         g1 = 1.d0
         g2 = 1.d0
         g3 = 1.d0
         g4 = 1.d0
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
         f1 = 1
         f2 = 3
         f3 = 2
         f4 = 4
         coef(1) = 1.d0
         coef(2) = -1.d0
         coef(3) = 1.d0
         coef(4) = 1.d0
         g1 = 1.d0
         g2 = 1.d0
         g3 = -1.d0
         g4 = 1.d0
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
      endif
!cccccccccccc slope limiter  cccccccccccccccc
    !     do j=1-cy, ncell_loc(2)+cy
    !       do i=1-cx, ncell_loc(1)+cx
         do j=0, ncell_loc(2)+1
           do i=0, ncell_loc(1)+1
           rhoL = g1*q_loc(i-cx,j-cy,f1)
           rhoM = g1*q_loc(i   ,j   ,f1)
           rhoR = g1*q_loc(i+cx,j+cy,f1)

            vxL = g2*q_loc(i-cx,j-cy,f2)/rhoL
            vxM = g2*q_loc(i   ,j   ,f2)/rhoM
            vxR = g2*q_loc(i+cx,j+cy,f2)/rhoR

            vyL = g3*q_loc(i-cx,j-cy,f3)/rhoL
            vyM = g3*q_loc(i   ,j   ,f3)/rhoM
            vyR = g3*q_loc(i+cx,j+cy,f3)/rhoR

           eneL = g4*q_loc(i-cx,j-cy,f4)/rhoL
           eneM = g4*q_loc(i   ,j   ,f4)/rhoM
           eneR = g4*q_loc(i+cx,j+cy,f4)/rhoR
  
             pL = (gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2))
             pM = (gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2))
             pR = (gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2))
#ifdef MHD
            vzL = g4*q_loc(i-cx,j-cy,f4)/rhoL
            vzM = g4*q_loc(i   ,j   ,f4)/rhoM
            vzR = g4*q_loc(i+cx,j+cy,f4)/rhoR

            bxL = 0.5d0*(g5*q_loc(i-cx,j-cy,f5)+g9*q_loc(i-cx,j-cy,f9))
            bxM = 0.5d0*(g5*q_loc(i   ,j   ,f5)+g9*q_loc(i   ,j   ,f9))
            bxR = 0.5d0*(g5*q_loc(i+cx,j+cy,f5)+g9*q_loc(i+cx,j+cy,f9))
 
            byL = 0.5d0*(g6*q_loc(i-cx,j-cy,f6)+g10*q_loc(i-cx,j-cy,f10))
            byM = 0.5d0*(g6*q_loc(i   ,j   ,f6)+g10*q_loc(i   ,j   ,f10))
            byR = 0.5d0*(g6*q_loc(i+cx,j+cy,f6)+g10*q_loc(i+cx,j+cy,f10))

            bzL = 0.5d0*(g7*q_loc(i-cx,j-cy,f7)+g11*q_loc(i-cx,j-cy,f11))
            bzM = 0.5d0*(g7*q_loc(i   ,j   ,f7)+g11*q_loc(i   ,j   ,f11))
            bzR = 0.5d0*(g7*q_loc(i+cx,j+cy,f7)+g11*q_loc(i+cx,j+cy,f11))

           eneL = g8*q_loc(i-cx,j-cy,f8)/rhoL
           eneM = g8*q_loc(i   ,j   ,f8)/rhoM
           eneR = g8*q_loc(i+cx,j+cy,f8)/rhoR

             pL = (gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2)-0.5d0*(bxL**2+byL**2+bzL**2))
             pM = (gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2+vzM**2)-0.5d0*(bxM**2+byM**2+bzM**2))
             pR = (gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2)-0.5d0*(bxR**2+byR**2+bzR**2))
#endif


         SL(i,j,1) = fslop(rhoL,rhoM,rhoR)
         SL(i,j,2) = fslop( vxL, vxM, vxR)
         SL(i,j,3) = fslop( vyL, vyM, vyR)
         SL(i,j,4) = fslop(  pL,  pM,  pR)
#ifdef MHD
         SL(i,j,4) = fslop( vzL, vzM, vzR)
         SL(i,j,5) = fslop( bxL, bxM, bxR)
         SL(i,j,6) = fslop( byL, byM, byR)
         SL(i,j,7) = fslop( bzL, bzM, bzR)
         SL(i,j,8) = fslop(  pL,  pM,  pR)
#endif
           enddo
         enddo

!cccccccccccc flux calculation cccccccccccccc
   !      do j=1-cy, ncell_loc(2)
   !        do i=1-cx, ncell_loc(1)
         do j=0, ncell_loc(2)+1
           do i=0, ncell_loc(1)+1

                   ql(1) = g1*q_loc(i,j,f1)
                   ql(2) = g2*q_loc(i,j,f2)/ql(1)
                   ql(3) = g3*q_loc(i,j,f3)/ql(1)
                   ql(4) = g4*(gam-1.d0)*(q_loc(i,j,f4)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2))

                   qr(1) = g1*q_loc(i+cx,j+cy,f1)
                   qr(2) = g2*q_loc(i+cx,j+cy,f2)/qr(1)
                   qr(3) = g3*q_loc(i+cx,j+cy,f3)/qr(1)
                   qr(4) = g4*(gam-1.d0)*(q_loc(i+cx,j+cy,f4)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2))

#ifdef MHD
                   ql(4) = g4*q_loc(i,j,f4)/ql(1)   ! vz
                   ql(5) = 0.5d0*(g5*q_loc(i,j,f5)+ g9*q_loc(i,j,f9 ))  ! bx
                   ql(6) = 0.5d0*(g6*q_loc(i,j,f6)+g10*q_loc(i,j,f10))  ! by
                   ql(7) = 0.5d0*(g7*q_loc(i,j,f7)+g11*q_loc(i,j,f11))  ! bz
                   ql(8) = g8*(gam-1.d0)*(q_loc(i,j,f8)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2+ql(4)**2)-0.5d0*(ql(5)**2+ql(6)**2+ql(7)**2))  !pressure

                   qr(4) = g4*q_loc(i+cx,j+cy,f4)/qr(1)   ! vz
                   qr(5) = 0.5d0*(g5*q_loc(i+cx,j+cy,f5)+ g9*q_loc(i+cx,j+cy,f9 )) !bx
                   qr(6) = 0.5d0*(g6*q_loc(i+cx,j+cy,f6)+g10*q_loc(i+cx,j+cy,f10)) !by
                   qr(7) = 0.5d0*(g7*q_loc(i+cx,j+cy,f7)+g11*q_loc(i+cx,j+cy,f11)) !bz
                   qr(8) = g8*(gam-1.d0)*(q_loc(i+cx,j+cy,f8)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2+qr(4)**2)-0.5d0*(qr(5)**2+qr(6)**2+qr(7)**2)) !pressure
#endif

               slopeL(:) = SL(i,j,:)
               slopeR(:) = SL(i+cx,j+cy,:)
!!!!!! ccccccc puting if-endif into loop is very bad
!!!!!! try to use preprocessor instead !!!!
!           if(solver_type .eq. 1) then
!              call fluxRieIso2d(ql,qr,slopeL,slopeR,Ftemp)
!           elseif(solver_type .eq. 2) then
!              call fluxHLLC2d(ql,qr,slopeL,slopeR,Ftemp)
!           endif
#ifdef ISOTHERMAL
#ifdef ExactRie
  call fluxRieIso2d(ql,qr,slopeL,slopeR,Ftemp)
#endif
#ifdef HLL
  call fluxHLLIso2d(ql,qr,slopeL,slopeR,Ftemp)
#endif
#endif

#ifdef ADIABATIC
#ifdef HLLC
  call fluxHLLC2d(ql,qr,slopeL,slopeR,Ftemp)
#endif
#ifdef HLL
#ifdef MHD
  call fluxHLLAdiMHD(ql,qr,slopeL,slopeR,Ftemp)
#else
  call fluxHLLAdi2d(ql,qr,slopeL,slopeR,Ftemp)
#endif
#endif
#ifdef HLLD
  call fluxHLLDAdiMHD(ql,qr,slopeL,slopeR,Ftemp)
#endif
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 Fx(i,j,:) = Ftemp(:)
          enddo
        enddo

       call bndflux2d(q_loc,Fx,dd)  ! implement your boundary flux in this subroutine

        dtddx = dt/dx
          do j=1, ncell_loc(2)
            do i=1, ncell_loc(1)
               temp2_loc(i,j,1) = temp1_loc(i,j,1)+ &
                 coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
               temp2_loc(i,j,2) = temp1_loc(i,j,2)+ &
                 coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
               temp2_loc(i,j,3) = temp1_loc(i,j,3)+ &
                 coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
               temp2_loc(i,j,4) = temp1_loc(i,j,4)+ &
                 coef(4)*(Fx(i-cx,j-cy,f4)-Fx(i,j,f4))*dtddx
#ifdef MHD
               !temp2_loc(i,j,5) = temp1_loc(i,j,5)+ &
               !  coef(5)*(Fx(i-cx,j-cy,f5)-Fx(i,j,f5))*dtddx
               !temp2_loc(i,j,6) = temp1_loc(i,j,6)+ &
               !  coef(6)*(Fx(i-cx,j-cy,f6)-Fx(i,j,f6))*dtddx
               temp2_loc(i,j,7) = temp1_loc(i,j,7)+ &
                 coef(7)*(Fx(i-cx,j-cy,f7)-Fx(i,j,f7))*dtddx
               temp2_loc(i,j,8) = temp1_loc(i,j,8)+ &
                 coef(8)*(Fx(i-cx,j-cy,f8)-Fx(i,j,f8))*dtddx

               !temp2_loc(i,j,9)  = temp2_loc(i,j,5)
               !temp2_loc(i,j,10) = temp2_loc(i,j,6)
               temp2_loc(i,j,11) = temp2_loc(i,j,7)
#endif
             enddo
          enddo

!!!!! Evolve B field
#ifdef MHD
      allocate(Ex(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))
      allocate(Ey(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))
      allocate(Ez(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))

      Ex = 0.d0
      Ey = 0.d0
      Ez = 0.d0
      

         do j=0,ncell_loc(2)
           do i=0,ncell_loc(1)
             if(dd .eq. 1) then
               Ez(i,j)=0.25d0*(-Fx(i,j,f6)-Fx(i,j+1,f6))*coef(6)
               Ey(i,j)=0.25d0*( Fx(i,j,f7)+Fx(i,j+1,f7))*coef(7)
              endif
             if(dd .eq. 2) then
               Ez(i,j)=0.25d0*(Fx(i+1,j,f5)+Fx(i,j,f5))*coef(5)
               Ex(i,j)=0.25d0*(-Fx(i+1,j,f7)-Fx(i,j,f7))*coef(7)
             endif
           enddo
         enddo
  do j=1,ncell_loc(2)
    do i=1,ncell_loc(1)
       temp2_loc(i,j,5)  = temp1_loc(i,j,5)+ &
            (Ez(i-1,j-1)-Ez(i-1,j))*dtddx
       temp2_loc(i,j,9)  = temp1_loc(i,j,9)+ &
            (Ez(i  ,j-1)-Ez(i  ,j))*dtddx


       temp2_loc(i,j,6)=temp1_loc(i,j,6)+ &
            (Ez(i,j-1)-Ez(i-1,j-1))*dtddx
       temp2_loc(i,j,10) = temp1_loc(i,j,10)+ &
            (Ez(i,j  )-Ez(i-1,j  ))*dtddx

       !temp2_loc(i,j,7)=temp1_loc(i,j,7)+ &
       !     (Ey(i-1,j-1)-Ey(i,j-1))*dtddx+ &
       !     (Ex(i-1,j)-Ex(i-1,j-1))*dtddx
       !temp2_loc(i,j,11) = temp2_loc(i,j,7)
            

    enddo
  enddo

deallocate(Ex)
deallocate(Ey)
deallocate(Ez)
#endif
deallocate(SL)
deallocate(Fx)
end subroutine
