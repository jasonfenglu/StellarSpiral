module plotting
REAL,SAVE                               ::points(4,2)
CONTAINS

SUBROUTINE StellarGasPlot(F,n,domain,contourn)
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)             ::F(:,:,:)        !plotting data
DOUBLE PRECISION,INTENT(IN)             ::domain        !plot range
INTEGER,INTENT(IN)                      ::contourn      !countour number
INTEGER,INTENT(IN)                      ::n             !gridsize
INTEGER                                 ::PGBEG
!open plotting device
IF (PGBEG(0,'/xs',1,1) .NE. 1)THEN
        write(6,*)achar(27)//'[33m pgplot open failed',achar(27)//'[0m'
        STOP
ENDIF
!call plotting routine
CALL output
IF (PGBEG(0,'ProjectedGasDensity.png/png',1,1) .NE. 1)THEN
        write(6,*)achar(27)//'[33m pgplot open failed',achar(27)//'[0m'
        STOP
ENDIF
!call plotting routine
CALL output

CONTAINS

SUBROUTINE output
IMPLICIT NONE
REAL                            ::dx,dy
REAL                            ::TR(6)
REAL                            ::BRIGHT,CONTRA
REAL                            ::vmax,vmin
REAL                            ::ALEV
INTEGER                         ::I
CALL PGSVP(0.0,0.95,0.0,0.95)
dx = real(domain)/real(n/2)
dy = real(domain)/real(n/2)
TR(3) = 0.
TR(5) = 0.
TR(2) = dx
TR(1) = -domain-dx/2.d0
TR(4) = -domain-dy/2.d0
TR(6) = dy

BRIGHT = 0.5
CONTRA = 0.9

vmax = real(MAXVAL(F(:,:,2)))
vmin = real(MINVAL(F(:,:,2)))
vmax = vmax * 1.1d0
vmin = vmin * 1.1d0
write(6,*)achar(27)//'[33m Ploting z scale :',vmax,vmin,achar(27)//'[0m'

!setting color style
CALL PALETT(2,CONTRA,Bright)
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)

!CALL PGSCI(3)
!DO I = 1, contourn
!        ALEV = vmin + (I-1)*(vmax - vmin)/real(contourn)
!        print *,I,ALEV
!        CALL PGCONS(real(F),N,N,1,N,1,N,ALEV,-1,TR)
!ENDDO
vmax = real(MAXVAL(F(:,:,2)))
vmin = real(MINVAL(F(:,:,2)))
vmax = 30.0
CALL PGIMAG(REAL(F(:,:,2)),n,n,1,n,1,n,vmin,vmax,TR)
CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
CALL PGPT(1,(/0./),(/0./),2)
CALL PGSCI(1)
CALL PGBBUF
CALL PGCLOS


ENDSUBROUTINE
ENDSUBROUTINE

SUBROUTINE countour(F,n,domain,contourn)
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)             ::F(:,:)        !plotting data
DOUBLE PRECISION,INTENT(IN)             ::domain        !plot range
REAL                                    ::TR(6)         !plot geometry
REAL                                    ::dx,dy
REAL                                    ::BRIGHT,CONTRA
REAL                                    ::ALEV(1)       !drawing elevation
REAL                                    ::vmax,vmin     !range of z-direction
INTEGER,INTENT(IN)                      ::n             !dimentsion
INTEGER                                 ::PGBEG
INTEGER                                 ::contourn      !countour number
INTEGER                                 ::i             !iteration index

!open plotting device
IF (PGBEG(0,'/xs',1,1) .NE. 1)THEN
        write(6,*)achar(27)//'[33m pgplot open failed',achar(27)//'[0m'
        STOP
ENDIF
!call plotting routine
CALL output
IF (PGBEG(0,'countour.ps/vcps',1,1) .NE. 1)THEN
        write(6,*)achar(27)//'[33m pgplot open failed',achar(27)//'[0m'
        STOP
ENDIF
!call plotting routine
CALL output

CONTAINS

SUBROUTINE output
CALL PGSVP(0.0,0.95,0.0,0.95)
dx = real(domain)/real(n/2)
dy = real(domain)/real(n/2)
TR(3) = 0.
TR(5) = 0.
TR(2) = dx
TR(1) = -domain-dx/2.d0
TR(4) = -domain-dy/2.d0
TR(6) = dy

BRIGHT = 0.5
CONTRA = 0.9

vmax = real(MAXVAL(F(:,:)))
vmin = real(MINVAL(F(:,:)))
vmax = vmax * 1.1d0
vmin = vmin * 1.1d0
write(6,*)achar(27)//'[33m Ploting z scale :',vmax,vmin,achar(27)//'[0m'

!setting color style
CALL PALETT(2,CONTRA,Bright)
CALL PGBBUF
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
CALL PGBBUF

CALL PGSCI(3)
DO I = 1, contourn
        ALEV = vmin + (I-1)*(vmax - vmin)/real(contourn)
        print *,I,ALEV
        CALL PGCONS(real(F),N,N,1,N,1,N,ALEV,-1,TR)
ENDDO
!CALL PGIMAG(REAL(F(:,:)),n,n,1,n,1,n,vmin,vmax,TR)
CALL PGPT(1,(/0./),(/0./),2)
CALL PGSCI(1)
CALL PGBBUF
CALL PGCLOS


ENDSUBROUTINE

ENDSUBROUTINE

SUBROUTINE plotdensity(F,n,domain)
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)             ::F(:,:)        !plotting data
DOUBLE PRECISION                        ::domain        !plot range
REAL                                    ::TR(6)         !plot geometry
REAL                                    ::TR2(6)        !plot geometry
REAL                                    ::vmax,vmin
REAL                                    ::BRIGHT,CONTRA
INTEGER                                 ::m,n           !dimentsion
INTEGER                                 ::PGBEG
REAL                                    ::dx,dy
INTEGER                                 ::noutput


IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
CALL output
IF (PGBEG(0,'density.png/png',1,1) .NE. 1) STOP
CALL output

CONTAINS

SUBROUTINE output()
LOGICAL                                 ::rauto,drawcir
REAL                                    ::den
namelist /plotpara1/ rauto,den,drawcir

open(20,file='para.list')
read(20,nml=plotpara1)
close(20)

CALL PGSVP(0.0,0.95,0.0,0.95)
m = n
dx = real(domain)/real(n)*2.
dy = real(domain)/real(m)*2.

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
vmax = vmax * 1.1d0
vmin = vmin * 1.1d0
write(6,*)achar(27)//'[33m Ploting z scale :',vmax,vmin,achar(27)//'[0m'
if(.not.rauto)then
        vmax = den
        vmin = -vmax
        write(6,*)achar(27)//'[33m Changed to z scale :',vmax,vmin,achar(27)//'[0m'
endif

CALL PALETT(2,CONTRA,Bright)
CALL PGBBUF
CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
CALL PGIMAG(REAL(F(:,:)),m,n,1,n,1,m,vmin,vmax,TR)
CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
CALL PGSCH(1.0)
CALL PGLAB('kpc','kpc','Density')
CALL PGSFS(2)
CALL PGSCI(0)
!CALL PGPT(4,points(:,1),points(:,2),2)
if(drawcir)then
        CALL PGCIRC(0.,0.,1.26)
        CALL PGCIRC(0.,0.,2.36)
        CALL PGCIRC(0.,0.,4.72)
        CALL PGCIRC(0.,0.,10.636)
        CALL PGCIRC(0.,0.,8.83)
endif

!CALL PGLINE(2,(/2.,7./),(/-8.,-8/))

CALL PGCLOS

ENDSUBROUTINE

ENDSUBROUTINE

        SUBROUTINE plotdensity2(F,n,domain,r,k3,kn,u)
        IMPLICIT NONE
        DOUBLE PRECISION,ALLOCATABLE,INTENT(IN) ::F(:,:,:)!plotting data
        DOUBLE PRECISION                        ::domain!plot range
        DOUBLE PRECISION                        ::r(:),k3(:,:),u(:,:)
        REAL                                    ::TR(6) !plot geometry
        REAL                                    ::TR2(6) !plot geometry
        REAL                                    ::vmax,vmin
        REAL                                    ::BRIGHT,CONTRA
        INTEGER                                 ::m,n,kn!dimentsion
        INTEGER                                 ::PGBEG
        REAL                                    ::dx,dy
        INTEGER                                 ::noutput


        IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
        CALL output
        IF (PGBEG(0,'density.png/png',1,1) .NE. 1) STOP
        CALL PGPAP(23.,0.618)
        CALL output
        IF (PGBEG(0,'2density.png/png',1,1) .NE. 1) STOP
        CALL PGPAP(20.,0.618)
        CALL outputd


        
        CONTAINS
        SUBROUTINE output()
        LOGICAL                                 ::rauto(2),drawcir(2)
        REAL                                    ::den(2),alpha
        namelist /plotpara/ rauto,den,drawcir,alpha

        open(10,file='para.list')
        read(10,nml=plotpara)
        close(10)

        CALL PGSVP(0.0,0.95,0.0,0.95)
        CALL PGSUBP(-3,2)
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

        CALL PALETT(2,CONTRA,Bright)
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        !!Density
        vmax = real(MAXVAL(F(:,:,1)))
        vmin = real(MINVAL(F(:,:,1)))
        vmax = vmax * 1.1d0
        vmin = vmin * 1.1d0
        print *,vmax,vmin
        if(.not.rauto(1))then
                        vmax = den(1)
                        vmin = -vmax
        endif
        CALL PGIMAG(REAL(F(:,:,1)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        if(drawcir(1))then
                CALL PGSFS(2)
                CALL PGSCI(0)
                CALL PGCIRC(0.,0.,1.26)
                CALL PGCIRC(0.,0.,2.36)
                CALL PGCIRC(0.,0.,4.72)
                CALL PGCIRC(0.,0.,10.636)
        endif
        CALL PGSCI(1)
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        !!Density
        vmax = real(MAXVAL(F(:,:,2)))
        vmin = real(MINVAL(F(:,:,2)))
        vmax = vmax * 1.1d0
        vmin = vmin * 1.1d0
        print *,vmax,vmin
        if(.not.rauto(2))then
                        vmax = den(2)
                        vmin = -vmax
        endif
        CALL PGIMAG(REAL(F(:,:,2)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        CALL PGSCH(1.0)
        CALL PGLAB('kpc','kpc','Density')
        CALL PGSFS(2)
        CALL PGSCI(0)
        CALL PGPT(4,points(:,1),points(:,2),2)
        if(drawcir(2))then
                CALL PGCIRC(0.,0.,1.26)
                CALL PGCIRC(0.,0.,2.36)
                CALL PGCIRC(0.,0.,4.72)
                CALL PGCIRC(0.,0.,10.636)
        endif

!       CALL PGLINE(2,(/2.,7./),(/-8.,-8/))

        !!k3
        CALL PGSCI(1)
        CALL PGENV(0.,real(domain),-1.,2.,0,1)
        CALL PGLINE(kn,real(r(:)),real(k3(:,1)))
        CALL PGSCI(2)
        CALL PGLINE(kn,real(r(:)),real(k3(:,2)))

        !!u
        CALL PGSCI(1)
        CALL PGENV(0.,real(domain),0.,8.,0,1)
        CALL PGLINE(kn,real(r(:)),real(u(:,1)))
        CALL PGSCI(2)
        CALL PGLINE(kn,real(r(:)),real(u(:,2)))
!       CALL PGLINE(kn,real(r(:)),real(d(:,1)))
!       CALL PGLINE(kn,real(r(:)),real(d(:,2)))


        !!map
        CALL PGSCI(1)
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        !!Density
        vmax = real(MAXVAL(F(:,:,1))+MAXVAL(F(:,:,2)))
        vmin = -vmax
        vmax = vmax * 1.1d0
        vmin = vmin * 1.1d0
        print *,vmax,vmin
        if(.not.rauto(1))then
                        vmax = den(1)
                        vmin = -vmax
        endif
        CALL PGIMAG(REAL(F(:,:,1))+REAL(F(:,:,2))*alpha,2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        if(drawcir(1))then
                CALL PGSFS(2)
                CALL PGSCI(0)
                CALL PGCIRC(0.,0.,1.26)
                CALL PGCIRC(0.,0.,2.36)
                CALL PGCIRC(0.,0.,4.72)
                CALL PGCIRC(0.,0.,10.636)
        endif
        CALL PGSCI(0)
        CALL PGPT(4,points(:,1),points(:,2),2)
        CALL PGSCI(1)

        !summantion of u
        CALL PGENV(0.,real(domain),0.,8.,0,1)
        CALL PGLINE(kn,real(r(:)),real(u(:,1)+u(:,2)*alpha))

        CALL PGCLOS

        ENDSUBROUTINE

        SUBROUTINE outputd()
        LOGICAL                                 ::rauto(2),drawcir(2)
        REAL                                    ::den(2),alpha
        namelist /plotpara/ rauto,den,drawcir,alpha

        open(10,file='para.list')
        read(10,nml=plotpara)
        close(10)

        CALL PGSVP(0.0,0.95,0.0,0.95)

        !!map
        CALL PALETT(2,CONTRA,Bright)
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        !!Density
        vmax = real(MAXVAL(F(:,:,1))+MAXVAL(F(:,:,2)))
        vmin = -vmax
        vmax = vmax * 1.1d0
        vmin = vmin * 1.1d0
        print *,vmax,vmin
        if(.not.rauto(1))then
                        vmax = den(1)
                        vmin = -vmax
        endif
        CALL PGIMAG(REAL(F(:,:,2)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        if(drawcir(1))then
                CALL PGSFS(2)
                CALL PGSCI(0)
                CALL PGCIRC(0.,0.,1.26)
                CALL PGCIRC(0.,0.,2.36)
                CALL PGCIRC(0.,0.,4.72)
                CALL PGCIRC(0.,0.,10.636)
        endif
        CALL PGSCI(0)
!       CALL PGPT(4,points(:,1),points(:,2),2)
        CALL PGSCI(1)

        CALL PGCLOS

        ENDSUBROUTINE
        ENDSUBROUTINE

        SUBROUTINE plotdensity3(F,n,domain,r,k3,kn,u)
        IMPLICIT NONE
        DOUBLE PRECISION,ALLOCATABLE,INTENT(IN) ::F(:,:,:)!plotting data
        DOUBLE PRECISION                        ::domain!plot range
        DOUBLE PRECISION                        ::r(:),k3(:,:),u(:,:)
        REAL                                    ::TR(6) !plot geometry
        REAL                                    ::TR2(6) !plot geometry
        REAL                                    ::vmax,vmin
        REAL                                    ::BRIGHT,CONTRA
        INTEGER                                 ::m,n,kn!dimentsion
        INTEGER                                 ::PGBEG
        REAL                                    ::dx,dy
        INTEGER                                 ::noutput


        IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
        CALL output
        IF (PGBEG(0,'density.png/png',1,1) .NE. 1) STOP
        CALL PGPAP(23.,0.618)
        CALL output
        IF (PGBEG(0,'2density.png/png',1,1) .NE. 1) STOP
        CALL PGPAP(20.,0.618)
        CALL outputd

        
        CONTAINS
        SUBROUTINE output()
        LOGICAL                                 ::rauto(2),drawcir(2)
        REAL                                    ::den(2),alpha
        namelist /plotpara/ rauto,den,drawcir,alpha

        open(10,file='para.list')
        read(10,nml=plotpara)
        close(10)

        CALL PGSVP(0.0,0.95,0.0,0.95)
        CALL PGSUBP(-3,2)
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

        CALL PALETT(2,CONTRA,Bright)
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        !!Density
        vmax = real(MAXVAL(F(:,:,1)))
        vmin = real(MINVAL(F(:,:,1)))
        vmax = vmax * 1.1d0
        vmin = vmin * 1.1d0
        print *,vmax,vmin
        if(.not.rauto(1))then
                        vmax = den(1)
                        vmin = -vmax
        endif
        CALL PGIMAG(REAL(F(:,:,1)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        if(drawcir(1))then
                CALL PGSFS(2)
                CALL PGSCI(0)
                CALL PGCIRC(0.,0.,1.26)
                CALL PGCIRC(0.,0.,2.36)
                CALL PGCIRC(0.,0.,4.72)
                CALL PGCIRC(0.,0.,10.636)
        endif
        CALL PGSCI(1)
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        !!Density
        vmax = real(MAXVAL(F(:,:,2)))
        vmin = real(MINVAL(F(:,:,2)))
        vmax = vmax * 1.1d0
        vmin = vmin * 1.1d0
        print *,vmax,vmin
        if(.not.rauto(2))then
                        vmax = den(2)
                        vmin = -vmax
        endif
        CALL PGIMAG(REAL(F(:,:,2)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        CALL PGSCH(1.0)
        CALL PGLAB('kpc','kpc','Density')
        CALL PGSFS(2)
        CALL PGSCI(0)
        CALL PGPT(4,points(:,1),points(:,2),2)
        if(drawcir(2))then
                CALL PGCIRC(0.,0.,1.26)
                CALL PGCIRC(0.,0.,2.36)
                CALL PGCIRC(0.,0.,4.72)
                CALL PGCIRC(0.,0.,10.636)
        endif

!       CALL PGLINE(2,(/2.,7./),(/-8.,-8/))

        !!k3
        CALL PGSCI(1)
        CALL PGENV(0.,real(domain),-1.,2.,0,1)
        CALL PGLINE(kn,real(r(:)),real(k3(:,1)))
        CALL PGSCI(2)
        CALL PGLINE(kn,real(r(:)),real(k3(:,2)))
        CALL PGSCI(3)
        CALL PGLINE(kn,real(r(:)),real(k3(:,3)))

        !!u
        CALL PGSCI(1)
        CALL PGENV(0.,real(domain),0.,8.,0,1)
        CALL PGLINE(kn,real(r(:)),real(u(:,1)))
        CALL PGSCI(2)
        CALL PGLINE(kn,real(r(:)),real(u(:,2)))
        CALL PGSCI(3)
        CALL PGLINE(kn,real(r(:)),real(u(:,3)))


        !!map
        CALL PGSCI(1)
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        !!Density
        vmax = real(MAXVAL(F(:,:,1))+MAXVAL(F(:,:,2)))
        vmin = -vmax
        vmax = vmax * 1.1d0
        vmin = vmin * 1.1d0
        print *,vmax,vmin
        if(.not.rauto(1))then
                        vmax = den(1)
                        vmin = -vmax
        endif
        CALL PGIMAG(REAL(F(:,:,1))+REAL(F(:,:,2))*alpha,2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        if(drawcir(1))then
                CALL PGSFS(2)
                CALL PGSCI(0)
                CALL PGCIRC(0.,0.,1.26)
                CALL PGCIRC(0.,0.,2.36)
                CALL PGCIRC(0.,0.,4.72)
                CALL PGCIRC(0.,0.,10.636)
        endif
        CALL PGSCI(0)
        CALL PGPT(4,points(:,1),points(:,2),2)
        CALL PGSCI(1)

        !summantion of u
        CALL PGENV(0.,real(domain),0.,8.,0,1)
        CALL PGLINE(kn,real(r(:)),real(u(:,1)+u(:,2)*alpha))

        CALL PGCLOS

        ENDSUBROUTINE

        SUBROUTINE outputd()
        LOGICAL                                 ::rauto(2),drawcir(2)
        REAL                                    ::den(2),alpha
        namelist /plotpara/ rauto,den,drawcir,alpha

        open(10,file='para.list')
        read(10,nml=plotpara)
        close(10)

        CALL PGSVP(0.0,0.95,0.0,0.95)

        !!map
        CALL PALETT(2,CONTRA,Bright)
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        !!Density
        vmax = real(MAXVAL(F(:,:,1))+MAXVAL(F(:,:,2)))
        vmin = -vmax
        vmax = vmax * 1.1d0
        vmin = vmin * 1.1d0
        print *,vmax,vmin
        if(.not.rauto(1))then
                        vmax = den(1)
                        vmin = -vmax
        endif
        CALL PGIMAG(REAL(F(:,:,1))+REAL(F(:,:,2))+REAL(F(:,:,3)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        if(drawcir(1))then
                CALL PGSFS(2)
                CALL PGSCI(0)
                CALL PGCIRC(0.,0.,1.26)
                CALL PGCIRC(0.,0.,2.36)
                CALL PGCIRC(0.,0.,4.72)
                CALL PGCIRC(0.,0.,10.636)
        endif
        CALL PGSCI(0)
        CALL PGPT(4,points(:,1),points(:,2),2)
        CALL PGSCI(1)

        CALL PGCLOS

        ENDSUBROUTINE
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
        CALL PGPAP(12.,0.618)

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
        vmax = 1.0
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
        
        SUBROUTINE plotpspdsearchset(F,n,m,domain)
        IMPLICIT NONE
        DOUBLE PRECISION                        ::F(:,:,:)!plotting data
        DOUBLE PRECISION                        ::domain(4)!plot range
        REAL                                    ::TR(6) !plot geometry
        REAL                                    ::TR2(6) !plot geometry
        REAL                                    ::vmax,vmin
        REAL                                    ::BRIGHT,CONTRA
        INTEGER                                 ::m,n   !dimentsion
        INTEGER                                 ::npanel,i
        INTEGER                                 ::PGBEG
        REAL                                    ::dx,dy

        IF (PGBEG(0,'searchallset.png/png',1,1) .NE. 1) STOP
        CALL PGPAP(23.,0.618)
!       IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
        CALL PGSVP(0.0,0.95,0.0,0.95)

        CALL PGSUBP(3,3)
        npanel = SIZE(F,3)


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


        vmax = 0.5
        vmin = 0.
        CALL PALETT(2,CONTRA,Bright)
        CALL PGBBUF
        DO i = 1, npanel
                CALL PGENV(real(domain(1)),real(domain(2)),real(domain(3)),real(domain(4)),0,0)
        CALL PGIMAG(REAL(F(:,:,i)),n,m,1,n,1,m,vmin,vmax,TR)
        ENDDO

!       CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        CALL PGSCH(1.0)
        CALL PGLAB('kpc','kpc','Absolute Value of Error')


        CALL PGCLOS
        ENDSUBROUTINE

        SUBROUTINE plotforce(F,F2,n,domain)
        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(IN)             ::F(:,:),F2(:,:)!plotting data
        DOUBLE PRECISION                        ::domain!plot range
        REAL                                    ::TR(6) !plot geometry
        REAL                                    ::TR2(6) !plot geometry
        REAL                                    ::vmax,vmin
        REAL                                    ::BRIGHT,CONTRA
        INTEGER                                 ::m,n   !dimentsion
        INTEGER                                 ::PGBEG
        REAL                                    ::dx,dy
        INTEGER                                 ::noutput


        IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
        CALL output
!       IF (PGBEG(0,'density.png/png',1,1) .NE. 1) STOP
!       CALL output
        
        CONTAINS
        SUBROUTINE output()
        LOGICAL                                 ::rauto,drawcir
        REAL                                    ::den
        namelist /plotpara1/ rauto,den,drawcir

        open(20,file='para.list')
        read(20,nml=plotpara1)
        close(20)

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
        vmax = vmax * 1.1d0
        vmin = vmin * 1.1d0
!       print *,vmax,vmin
!       if(.not.rauto)then
!                       vmax = den
!                       vmin = -vmax
!       endif
        CALL PALETT(2,CONTRA,Bright)
        CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
        CALL PGIMAG(REAL(F(:,:)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
        CALL PGSCH(1.0)
        CALL PGLAB('kpc','kpc','Density')
        CALL PGSFS(2)
        CALL PGSCI(0)
!       CALL PGPT(4,points(:,1),points(:,2),2)
!       if(drawcir)then
!               CALL PGCIRC(0.,0.,1.26)
!               CALL PGCIRC(0.,0.,2.36)
!               CALL PGCIRC(0.,0.,4.72)
!               CALL PGCIRC(0.,0.,10.636)
!       endif

!       CALL PGLINE(2,(/2.,7./),(/-8.,-8/))


        ENDSUBROUTINE

        ENDSUBROUTINE

SUBROUTINE meshplot(dat,m,domain,zmax,n_in,zmin_in)
IMPLICIT NONE
!!input arguments
DOUBLE PRECISION,INTENT(IN)             ::dat(:,:)
DOUBLE PRECISION,INTENT(IN)             ::domain
DOUBLE PRECISION,INTENT(IN)             ::zmax
DOUBLE PRECISION,INTENT(IN),OPTIONAL    ::zmin_in
DOUBLE PRECISION                        ::zmin
INTEGER,INTENT(IN)                      ::m
INTEGER,INTENT(IN),OPTIONAL             ::n_in
INTEGER                                 ::n
!!used variables
REAL                                    ::TR(6)         !plot geometry
REAL                                    ::dx,dy
REAL                                    ::BRIGHT,CONTRA

!check if zmin is passed 
if(present(zmin_in))then
        zmin = zmin_in
ELSE
        zmin = - zmax
ENDIF
!check if n is passed 
if(present(n_in))THEN
        n = n_in
ELSE
        n = m
ENDIF
dx = real(domain)/real(n/2)
dy = real(domain)/real(m/2)
TR(3) = 0.
TR(5) = 0.
TR(2) = dx
TR(1) = -real(domain)-dx/2.d0
TR(4) = -real(domain)-dy/2.d0
TR(6) = dy
BRIGHT = 0.5
CONTRA = 0.9
CALL PALETT(2,CONTRA,Bright)
CALL PGIMAG(REAL(dat(:,:)),m,n,1,n,1,m,real(zmin),real(zmax),TR)
CALL PGWEDG('RI', 1.0, 4.0, real(zmin), real(zmax), '')
ENDSUBROUTINE

SUBROUTINE contourplot(dat,m,domain,contourn,n_in)
IMPLICIT NONE
!!input arguments
DOUBLE PRECISION,INTENT(IN)             ::dat(:,:)
DOUBLE PRECISION,INTENT(IN)             ::domain
DOUBLE PRECISION                        ::zmax
DOUBLE PRECISION                        ::zmin_in
DOUBLE PRECISION                        ::zmin
REAL                                    ::ALEV(1)       !drawing elevation
INTEGER,INTENT(IN)                      ::m
INTEGER,INTENT(IN),OPTIONAL             ::n_in
INTEGER,INTENT(IN)                      ::contourn      !countour number
INTEGER                                 ::n
INTEGER                                 ::I
!!used variables
REAL                                    ::TR(6)         !plot geometry
REAL                                    ::dx,dy
REAL                                    ::BRIGHT,CONTRA

!check if n is passed 
if(present(n_in))THEN
        n = n_in
ELSE
        n = m
ENDIF
dx = real(domain)/real(n/2)
dy = real(domain)/real(m/2)
TR(3) = 0.
TR(5) = 0.
TR(2) = dx
TR(1) = -real(domain)-dx/2.d0
TR(4) = -real(domain)-dy/2.d0
TR(6) = dy
BRIGHT = 0.5
CONTRA = 0.9

zmax = maxval(dat)
zmin = minval(dat)

CALL PALETT(2,CONTRA,Bright)
DO I = 1, contourn
        ALEV = zmin + (I-1)*(zmax - zmin)/real(contourn)
        print *,I,ALEV
        CALL PGCONS(real(dat),N,N,1,N,1,N,ALEV,-1,TR)
ENDDO
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
