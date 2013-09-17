PROGRAM StellarGas
USE STELLARDISK
USE io
IMPLICIT NONE
CHARACTER(len=40)                       ::densityh5
CHARACTER(len=40)                       ::arg
type(typspiral)                         ::spiral
DOUBLE PRECISION,ALLOCATABLE            ::f(:,:,:)
DOUBLE PRECISION,ALLOCATABLE            ::x(:),y(:)
DOUBLE PRECISION,ALLOCATABLE            ::lf(:,:)
DOUBLE PRECISION                        ::r,th
INTEGER                                 ::fsize(2)
INTEGER                                 ::i,j

if(iargc().eq.1)then
        CALL getarg(1,arg)
        READ(arg,*)densityh5
        print *,'reading file: ',densityh5
        else 
                write(0,*)'usage: stellar filename.h5'
                STOP
endif

CALL readdensity


CALL stdpara.readstd
CALL spiral.init(500,12.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral)

!$OMP PARALLEL SHARED(f,spiral) PRIVATE(j,r,th)
!$OMP DO 
DO i = 1, fsize(1)
DO j = 1, fsize(1)
        r = sqrt(x(i)**2+y(j)**2)
        th = atan2(y(j),x(i))
        f(i,j,4) = max(sigma1(r,th,spiral),0.d0)
ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL 


CALL plotdensity

CONTAINS
SUBROUTINE readdensity
INTEGER,PARAMETER                       ::buf = 2
INTEGER                                 ::ndim(2)
ndim = h5size(densityh5,"density",2)

ALLOCATE(f(ndim(1),ndim(2),4))
CALL h5read(f(:,:,1),ndim(1),ndim(2),densityh5,"density")
CALL h5read(f(:,:,2),ndim(1),ndim(2),densityh5,"momx")
CALL h5read(f(:,:,3),ndim(1),ndim(2),densityh5,"momy")

ALLOCATE(x(ndim(1)))
ALLOCATE(y(ndim(1)))
CALL h5read(x,size(x),densityh5,"x")
CALL h5read(y,size(x),densityh5,"y")

fsize = ndim

ENDSUBROUTINE

SUBROUTINE plotdensity
IMPLICIT NONE
INTEGER                                 ::PGBEG

IF (PGBEG(0,'/xs',1,1) .NE. 1) STOP
CALL output
!IF (PGBEG(0,'density.png/png',1,1) .NE. 1) STOP
!CALL output

ENDSUBROUTINE

SUBROUTINE output()
IMPLICIT NONE
REAL                                    ::TR(6)         !plot geometry
REAL                                    ::vmax,vmin
REAL                                    ::BRIGHT,CONTRA
REAL                                    ::ALEV(1)       !drawing elevation
INTEGER,PARAMETER                       ::contourn=4    !countour number

CALL PGSVP(0.0,0.95,0.0,0.95)

TR(3) = 0.
TR(5) = 0.
TR(2) = real(-x(1) + x(2))
TR(1) = real(2.d0*x(1) - x(2))
TR(4) = real(2.d0*y(1) - y(2))
TR(6) = real(-y(1) + y(2))

BRIGHT = 0.5
CONTRA = 0.9


!!Density
vmax = real(MAXVAL(F(:,:,1)))
vmin = real(MINVAL(F(:,:,1)))
vmax = vmax * 1.1d0
vmin = vmin * 1.1d0
vmax = 3.
vmin = 0.
write(6,*)achar(27)//'[33m Ploting z scale :',vmax,vmin,achar(27)//'[0m'

CALL PALETT(2,CONTRA,Bright)
CALL PGBBUF
CALL PGENV(real(x(1)),-real(x(1)),real(y(1)),-real(y(1)),1,0)
CALL PGIMAG(REAL(F(:,:,1)),fsize(1),fsize(2),1,fsize(1),1,fsize(2),vmin,vmax,TR)
CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
CALL PGSCH(1.0)
CALL PGLAB('kpc','kpc','Density '//densityh5)
CALL PGSFS(2)
CALL PGSCI(0)


!!plot contour
vmax = real(MAXVAL(F(:,:,4)))
vmin = real(MINVAL(F(:,:,4)))
vmax = vmax * 1.1d0
vmin = vmin * 1.1d0
CALL PGSCI(13)
DO I = 1, contourn
        ALEV = vmin + (I-1)*(vmax - vmin)/real(contourn)
        CALL PGCONS(real(F(:,:,4)),fsize(1),fsize(1),1,fsize(1),1,fsize(1),ALEV,-1,TR)
ENDDO
CALL PGSCI(0)


!CALL PGPT(4,points(:,1),points(:,2),2)

!CALL PGLINE(2,(/2.,7./),(/-8.,-8/))

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

ENDPROGRAM
