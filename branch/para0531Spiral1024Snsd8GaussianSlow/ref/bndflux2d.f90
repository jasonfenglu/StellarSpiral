subroutine bndflux2d(q_loc,Fx,dd)
use common_params
implicit none
double precision:: q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,nvar)
double precision::    Fx(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,nvar)
double precision,dimension(:,:,:),allocatable:: Fbnd
integer         :: dd
integer         :: i


allocate(Fbnd(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,nvar))

if(myid .eq. 0 .and. dd .eq. 2) then
    call bndflux_outgoing_isothermal(q_loc, Fbnd, dd, snd)
         do i=0, ncell_loc(1)
            Fx(i,0,:) = Fbnd(i,0,:)
         enddo
endif

if(myid .eq. nprocs-1 .and. dd .eq. 2) then
    call bndflux_outgoing_isothermal(q_loc, Fbnd, dd, snd)
         do i=0,ncell_loc(1) 
            Fx(i,ncell_loc(2),:) = Fbnd(i,ncell_loc(2),:)
         enddo
endif


if(dd .eq. 1) then
    call bndflux_outgoing_isothermal(q_loc, Fx, dd)
endif

deallocate(Fbnd)
end subroutine
!!!==========================================================
       subroutine bndflux_outgoing_isothermal(U1, Fx, dd)
       use common_params
       implicit none
       dimension U1(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,nvar)
       dimension Fx(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,nvar)
       integer dd,i,j
       double precision:: csq,rho,ux,uy
       double precision:: r,z1,z2,z3,dW,dU,dV
       double precision:: U1,Fx,Wb,Ub,Vb
       double precision:: c, xl, yl,v0
       double precision:: pspd,A,B,maxrho
       pspd=21.6385d0
       maxrho=15.d0
       A=-0.249125d0
       B=0.552313d0

        dx = x_loc(2)-x_loc(1)
       csq = snd**2.d0
         c = snd
!       write(*,*) " snd:", snd
        v0 = 549.959d0
      if(dd .eq. 1) then
!cccccccccccccccccccccccccccccccccccccccccc
      do j = 1,ncell_loc(2)
         i   = 1
         r   = dsqrt(x_loc(i)**2.d0+y_loc(j)**2.d0)
         rho = maxrho*dexp(-(r/14.311d0)**2.d0)
         ! rho=1.53d0
          ux = -v0*y_loc(j)/(r**B+r**(1-A)) + y_loc(j)*pspd
          uy =  v0*x_loc(i)/(r**B+r**(1-A)) - x_loc(i)*pspd
          dW = U1(1,j,1)-rho
          dU = U1(1,j,2)-rho*ux
          dV = U1(1,j,3)-rho*uy
         !i   = 0
         xl = x_loc(i)-0.5d0*dx
         r   = dsqrt(xl**2.d0+y_loc(j)**2.d0)
          ux = -v0*y_loc(j)/(r**B+r**(1-A)) + y_loc(j)*pspd
          uy =  v0*xl/(r**B+r**(1-A)) - xl*pspd
          z1 =  ( ( ux+c )*dW-dU )/2.d0/c
          z2 =    (-uy   )*dW+dV
          z3 =  (-( ux-c )*dW+dU )/2.d0/c
         if ( (ux+c).gt.0.d0) z3 = 0.d0
         if ( (ux  ).gt.0.d0) z2 = 0.d0
         if ( (ux-c).gt.0.d0) z1 = 0.d0
         Wb      = rho   +(        z1+       z3    )
         Ub      = rho*ux+( (ux-c)*z1+(ux+c)*z3    )
         Vb      = rho*uy+( (uy  )*z1+(uy  )*z3+z2 )
         Fx(0,j,1) = Ub
         Fx(0,j,2) = Wb*( (Ub/Wb)**2.d0+csq )
         Fx(0,j,3) = Vb*(  Ub/Wb         )
      enddo

      do j = 1,ncell_loc(2)
         i   = ncell_loc(1)
         r   = dsqrt(x_loc(i)**2.d0+y_loc(j)**2.d0)
         rho = maxrho*dexp(-(r/14.311d0)**2.d0)
         !rho=1.53d0
          ux = -v0*y_loc(j)/(r**B+r**(1-A)) + y_loc(j)*pspd
          uy =  v0*x_loc(i)/(r**B+r**(1-A)) - x_loc(i)*pspd
          dW = U1(ncell_loc(1),j,1)-rho
          dU = U1(ncell_loc(1),j,2)-rho*ux
          dV = U1(ncell_loc(1),j,3)-rho*uy
          xl = x_loc(i)+0.5d0*dx
         r   = dsqrt(xl**2.d0+y_loc(j)**2.d0)
          ux = -v0*y_loc(j)/(r**B+r**(1-A)) + y_loc(j)*pspd
          uy =  v0*xl/(r**B+r**(1-A)) - xl*pspd
          z1 =  ( ( ux+c )*dW-dU )/2.d0/c
          z2 =    (-uy   )*dW+dV
          z3 =  (-( ux-c )*dW+dU )/2.d0/c
         if ( (ux+c).lt.0.d0) z3 = 0.d0
         if ( (ux  ).lt.0.d0) z2 = 0.d0
         if ( (ux-c).lt.0.d0) z1 = 0.d0
         Wb      = rho   +(        z1+       z3    )
         Ub      = rho*ux+( (ux-c)*z1+(ux+c)*z3    )
         Vb      = rho*uy+( (uy  )*z1+(uy  )*z3+z2 )
         Fx(i,j,1) = Ub
         Fx(i,j,2) = Wb*( (Ub/Wb)**2.d0+csq )
         Fx(i,j,3) = Vb*(  Ub/Wb         )
      enddo
!cccccccccccccccccccccccccccccccccccccccccccc
      elseif(dd .eq. 2) then
      do i = 1,ncell_loc(1)
         j   = 1
         r   = dsqrt(x_loc(i)**2.d0+y_loc(j)**2.d0)
         rho = maxrho*dexp(-(r/14.311d0)**2.d0)
         !rho=1.53d0
          uy = v0*y_loc(j)/(r**B+r**(1-A)) - y_loc(j)*pspd
          ux = v0*x_loc(i)/(r**B+r**(1-A)) - x_loc(i)*pspd
          dW = U1(i,1,1)-rho
          dU = U1(i,1,3)-rho*ux
          dV = -U1(i,1,2)-rho*uy
         !j   = 0
          yl = y_loc(j)-0.5d0*dx
         r   = dsqrt(x_loc(i)**2.d0+yl**2.d0)
          uy =  v0*yl/(r**B+r**(1-A)) - yl*pspd
          ux =  v0*x_loc(i)/(r**B+r**(1-A)) - x_loc(i)*pspd
          z1 =  ( ( ux+c )*dW-dU )/2.d0/c
          z2 =    (-uy   )*dW+dV
          z3 =  (-( ux-c )*dW+dU )/2.d0/c
         if ( (ux+c).gt.0.d0) z3 = 0.d0
         if ( (ux  ).gt.0.d0) z2 = 0.d0
         if ( (ux-c).gt.0.d0) z1 = 0.d0
         Wb      = rho   +(        z1+       z3    )
         Ub      = rho*ux+( (ux-c)*z1+(ux+c)*z3    )
         Vb      = rho*uy+( (uy  )*z1+(uy  )*z3+z2 )
         Fx(i,0,1) = Ub
         Fx(i,0,2) = Wb*( (Ub/Wb)**2+csq )
         Fx(i,0,3) = Vb*(  Ub/Wb         )
      enddo

      do i = 1,ncell_loc(1)
         j   = ncell_loc(2)
         r   = dsqrt(x_loc(i)**2.d0+y_loc(j)**2.d0)
         rho = maxrho*dexp(-(r/14.311d0)**2.d0)
         !rho=1.53d0
          uy = v0*y_loc(j)/(r**B+r**(1-A)) - y_loc(j)*pspd
          ux = v0*x_loc(i)/(r**B+r**(1-A)) - x_loc(i)*pspd
          dW = U1(i,ncell_loc(2),1)-rho
          dU = U1(i,ncell_loc(2),3)-rho*ux
          dV = -U1(i,ncell_loc(2),2)-rho*uy
          yl = y_loc(j)+0.5d0*dx
         r   = dsqrt(x_loc(i)**2.d0+yl**2.d0)
          uy =  v0*yl/(r**B+r**(1-A)) - yl*pspd
          ux =  v0*x_loc(i)/(r**B+r**(1-A)) - x_loc(i)*pspd
          z1 =  ( ( ux+c )*dW-dU )/2.d0/c
          z2 =    (-uy   )*dW+dV
          z3 =  (-( ux-c )*dW+dU )/2.d0/c
         if ( (ux+c).lt.0.d0) z3 = 0.d0
         if ( (ux  ).lt.0.d0) z2 = 0.d0
         if ( (ux-c).lt.0.d0) z1 = 0.d0
         Wb      = rho   +(        z1+       z3    )
         Ub      = rho*ux+( (ux-c)*z1+(ux+c)*z3    )
         Vb      = rho*uy+( (uy  )*z1+(uy  )*z3+z2 )
         Fx(i,j,1) = Ub
         Fx(i,j,2) = Wb*( (Ub/Wb)**2.d0+csq )
         Fx(i,j,3) = Vb*(  Ub/Wb         )
      enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        else
        write(*,*) " Ultra-dimension?!"
        endif
       end subroutine
