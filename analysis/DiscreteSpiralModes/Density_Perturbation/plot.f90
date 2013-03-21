        module plotting
        REAL,SAVE                               ::points(4,2)
        CONTAINS

        SUBROUTINE plotlog(dat,m,n)
        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(IN)             ::dat(:,:)
        DOUBLE PRECISION,ALLOCATABLE            ::plotrange(:,:)
        INTEGER                                 ::m,n   !dimentsion
        INTEGER                                 ::PGBEG,i
        
        
        IF (PGBEG(0,'./plotlog/log.png/png',1,1) .NE. 1) STOP
!       CALL PGSVP(0.0,0.95,0.0,0.95)
        
        
        ALLOCATE(plotrange(m,2))
        do i = 2, 4
                plotrange(i,1)=0.
                plotrange(i,2)=maxval(real(dat(i,:)))
        enddo
        plotrange(3,:) = (/-1.d0,5.d0/)

        do i = 2, 4
              CALL PGENV(0.,real(maxval(dat(1,:)))/2.,real(plotrange(i,1)),real(plotrange(i,2)),0,0)
              CALL PGLINE(n,real(dat(1,:)),real(dat(i,:)))
        enddo
        CALL PGCLOS
        ENDSUBROUTINE

        SUBROUTINE plotdensity(F,F2,force,n,domain)
        IMPLICIT NONE
        DOUBLE PRECISION,ALLOCATABLE,INTENT(IN) ::F(:,:),F2(:,:)!plotting data
        DOUBLE PRECISION                        ::force(:,:,:)
        DOUBLE PRECISION                        ::domain!plot range
        REAL                                    ::TR(6) !plot geometry
        REAL                                    ::TR2(6) !plot geometry
        REAL                                    ::vmax,vmin
        REAL                                    ::BRIGHT,CONTRA
        INTEGER                                 ::m,n   !dimentsion
        INTEGER                                 ::PGBEG
        REAL                                    ::dx,dy


        IF (PGBEG(0,'density.png/png',1,1) .NE. 1) STOP
!       IF (PGBEG(0,'/xserve',1,1) .NE. 1) STOP
        CALL PGSVP(0.0,0.95,0.0,0.95)

        m = n
        dx = real(domain)/real(n)
        dy = real(domain)/real(m)

        TR(3) = 0.
        TR(5) = 0.
        TR(2) = dx
        TR(1) = -domain-dx/2.d0
        TR(4) = -domain-dy/2.d0
        TR(6) = dy

        BRIGHT = 0.5
        CONTRA = 0.9


        !!Density
        vmax = real(MAXVAL(F(:,:)))
        vmin = real(MINVAL(F(:,:)))
        print *,vmax,vmin
!       vmax = vmax * 1.1d0
!       vmin = vmin * 1.1d0
        vmax = 350.
        vmin =-350.
        CALL PALETT(2,CONTRA,Bright)
        CALL PGBBUF
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        CALL PGIMAG(REAL(F(:,:)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        CALL PGSCH(1.0)
        CALL PGLAB('kpc','kpc','Density')
        CALL PGSFS(2)
        CALL PGSCI(0)
!       CALL PGPT(4,points(:,1),points(:,2),2)
!       CALL PGCIRC(0.,0.,1.26)
!       CALL PGCIRC(0.,0.,2.36)
!       CALL PGCIRC(0.,0.,4.72)
!       CALL PGCIRC(0.,0.,10.636)

        CALL PGLINE(2,(/2.,7./),(/-8.,-8/))

        !!Potential
!       vmax = real(MAXVAL(F2(:,:)))
!       vmin = real(MINVAL(F2(:,:)))
!       vmax = 100.
!       vmin = -100.
!       CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
!       CALL PGIMAG(REAL(F2),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
!       CALL PGWEDG('RI', 1.0, 4.0, vmax, vmin, '')
!       CALL PGSCH(1.0)
!       CALL PGLAB('kpc','kpc','Potential')

        !!Force
!       TR2 = 0.
!       TR2(2) = 8.d0*dx
!       TR2(1) = -domain-dx*4.d0
!       TR2(6) = 8.d0*dy
!       TR2(4) = -domain-dy*4.d0
!       CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,-1)
!       CALL PGIMAG(REAL(F2),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
!       CALL PGWEDG('RI', 1.0, 4.0, vmax, vmin, '')
!       CALL PGSCH(1.0)
!       CALL PGLAB('kpc','kpc','Force')
!       CALL PGSCH(0.8)
!       CALL PGSCI(0)
!       CALL PGSAH(1,20.,0.3)
!       CALL PGVECT(real(force(:,:,1)),real(force(:,:,2)),n/4,n/4,         &
!                   2,n/4-2, &
!                   2,n/4-2, &
!                   0.02,2,TR2,-1.E10)
        CALL PGCLOS
        ENDSUBROUTINE

        SUBROUTINE plotpspdsearch(F,n,m,domain)
        IMPLICIT NONE
        DOUBLE PRECISION                        ::F(:,:)!plotting data
        DOUBLE PRECISION                        ::domain(4)!plot range
        REAL                                    ::TR(6) !plot geometry
        REAL                                    ::TR2(6) !plot geometry
        REAL                                    ::vmax,vmin
        REAL                                    ::BRIGHT,CONTRA
        INTEGER                                 ::m,n   !dimentsion
        INTEGER                                 ::PGBEG
        REAL                                    ::dx,dy


        IF (PGBEG(0,'searchall.png/png',1,1) .NE. 1) STOP
        CALL PGSVP(0.0,0.95,0.0,0.95)

        dx = real(domain(2)-domain(1))/real(n)
        dy = real(domain(4)-domain(3))/real(m)

        TR(3) = 0.
        TR(5) = 0.
        TR(2) = REAL(domain(2)-domain(1))/REAL(n-1)
        TR(1) = REAL(domain(1))-TR(2)
        TR(6) = REAL(domain(4)-domain(3))/REAL(m-1)
        TR(4) = REAL(domain(3))-TR(6)

        BRIGHT = 0.5
        CONTRA = -0.9


!       vmax = real(MAXVAL(F(:,:)))
!       vmin = real(MINVAL(F(:,:)))
        vmax = 0.5
        vmin = 0.
        CALL PALETT(2,CONTRA,Bright)
        CALL PGBBUF
        CALL PGENV(real(domain(1)),real(domain(2)),real(domain(3)),real(domain(4)),0,0)
        CALL PGIMAG(REAL(F(:,:)),n,m,1,n,1,m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        CALL PGSCH(1.0)
        CALL PGLAB('kpc','kpc','Absolute Value of Error')


        CALL PGCLOS
        ENDSUBROUTINE

      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
!-----------------------------------------------------------------------
! Set a "palette" of colors in the range of color indices used by
! PGIMAG.
!-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
!
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
!
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
!
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
!
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
               0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
               0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
               0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
      IF (TYPE.EQ.1) THEN
!        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
!        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
!        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
!        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
!        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      ENDSUBROUTINE

        endmodule
