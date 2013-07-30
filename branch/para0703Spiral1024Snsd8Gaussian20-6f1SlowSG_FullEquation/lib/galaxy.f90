module galaxy
USE STELLARDISK
IMPLICIT NONE
type(typspiral),SAVE                   ::spiral
type    typspiralforce
        DOUBLE PRECISION,ALLOCATABLE   ::den(:,:)
        DOUBLE PRECISION,ALLOCATABLE   ::dden(:,:,:)
        double precision,dimension(:,:),allocatable::sgx,sgy !selfgravity
        double complex  ,dimension(:,:),allocatable::sgxker,sgyker !kernal 
endtype
type(typspiralforce),SAVE               ::stellarforce
DOUBLE PRECISION,ALLOCATABLE            ::InitGasDensity(:,:)
DOUBLE PRECISION,PARAMETER              ::GasDiskLength = 7.d0
DOUBLE PRECISION,PARAMETER              ::GasDiskSigma0 = 19.33d0

CONTAINS

FUNCTION RC(r)
IMPLICIT NONE
DOUBLE PRECISION        ::RC,r

RC = Omega(r,spiral)*r

ENDFUNCTION

FUNCTION gasdensity(r,z)
IMPLICIT NONE
DOUBLE PRECISION                ::gasdensity,r
DOUBLE PRECISION                ::z

gasdensity = GaussianDensity(r,z)
!gasdensity = UniDensity(r,z)
!gasdensity = Sigma0(r,spiral)
!gasdensity = Sigma1(r,th,spiral)

ENDFUNCTION

FUNCTION UniDensity(r,z)
IMPLICIT NONE
DOUBLE PRECISION        ::UniDensity,r,z

!Uniform density distribution
!
!
!       rho = const
!
!---------------------------

UniDensity = 1.d0

ENDFUNCTION
 
FUNCTION GaussianDensity(r,z)
IMPLICIT NONE
DOUBLE PRECISION                ::GaussianDensity,r
DOUBLE PRECISION,OPTIONAL       ::z
DOUBLE PRECISION                ::rho,a

!Gaussian density distribution
!
!
!  density = rho * exp(-r^2/(2a**2))
!
!---------------------------

rho = GasDiskSigma0
a   = GasDiskLength

GaussianDensity = rho * dexp(-r**2/2.d0/a**2)

ENDFUNCTION

FUNCTION dGaussianDensity(r,z)
IMPLICIT NONE
DOUBLE PRECISION        ::dGaussianDensity,r
DOUBLE PRECISION,OPTIONAL       ::z
DOUBLE PRECISION        ::rho,a

!Uniform density distribution
!
!
!  density = rho * exp(-r^2/(2a))
!
!---------------------------

rho = 1.d0
a   = 7.d0

dGaussianDensity = -1.d0*r/a*GaussianDensity(r,z)

ENDFUNCTION

subroutine write_stellargravity2d()
use common_params
implicit none
character(len=8)::flnm='force.h5'
character(len=20)::dsetname


!flnm( 1:1 ) = 'M'
!flnm( 2:2 ) = char(ichar('0')+ic1)
!flnm( 3:3 ) = char(ichar('0')+ic2)
!flnm( 4:4 ) = char(ichar('0')+ic3)
!flnm( 5:5 ) = char(ichar('0')+ic4)
!flnm( 6:8 ) = '.h5'

if(myid.eq.0)then
        write(6,*)achar(27)//'[95m writing spiral force'//achar(27)//'[0m'
ENDIF
dsetname='stellarfx'
call output2d(stellarforce.sgx,ncell(1),ncell(2),ncell_loc(1),ncell_loc(2),0,0,0,0,flnm,dsetname,1)
dsetname='stellarfy'
call output2d(stellarforce.sgy,ncell(1),ncell(2),ncell_loc(1),ncell_loc(2),0,0,0,0,flnm,dsetname,0)
dsetname='density'
call output2d(stellarforce.den,ncell(1),ncell(2),ncell_loc(1),ncell_loc(2),0,0,0,0,flnm,dsetname,0)
dsetname='x'
call output1d(x_loc,ncell(1),ncell_loc(1),ibuf,2,flnm,dsetname,0) !1: new 0:open
dsetname='y'
call output1d(y_loc,ncell(2),ncell_loc(2),jbuf,2,flnm,dsetname,0)
end subroutine

SUBROUTINE InitStellarDensity
USE common_params,only:ncell_loc,x_loc,y_loc
INTEGER                                 ::i,j
DOUBLE PRECISION                        ::x,y,r,th
DOUBLE PRECISION                        ::amp
namelist /FORCENML/ amp

open(10,file='para.list')
read(10,nml=FORCENML)
close(10)

DO j = 1, ncell_loc(2)
DO i = 1, ncell_loc(1)
        x = x_loc(i)
        y = y_loc(j)
        r = sqrt(x**2+y**2)
        th = atan2(y,x)
        stellarforce.den(i,j) = Sigma1(r,th,spiral)*amp
ENDDO
ENDDO
ENDSUBROUTINE

SUBROUTINE densityplot
!USE common_params,only:xrange
!DOUBLE PRECISION                ::domain
!INTEGER                         ::M,N
!
!domain = xrange(2)
!
!IF (PGBEG(0,'density.png/png',1,1) .NE. 1) STOP
!CALL PGSVP(0.0,0.95,0.0,0.95)
!CALL PGSUBP(-2,2)
!m = n
!dx = real(domain)/real(n)
!dy = real(domain)/real(m)
!
!TR(3) = 0.
!TR(5) = 0.
!TR(2) = dx
!TR(1) = -domain-dx/2.d0
!TR(4) = -domain-dy/2.d0
!TR(6) = dy
!
!BRIGHT = 0.5
!CONTRA = 0.9
!
!
!!!Density
!vmax = real(MAXVAL(F(:,:)))
!vmin = real(MINVAL(F(:,:)))
!vmax = vmax * 1.1d0
!vmin = vmin * 1.1d0
!print *,vmax,vmin
!if(.not.rauto)then
!        vmax = den
!        vmin = -vmax
!endif
!
!
!CALL PALETT(2,CONTRA,Bright)
!CALL PGBBUF
!CALL PGENV(-real(domain),real(domain),-real(domain),real(domain),1,0)
!CALL PGIMAG(REAL(F(:,:)),2*m,2*n,1,2*n,1,2*m,vmin,vmax,TR)
!CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, '')
!CALL PGSCH(1.0)
!
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
!      IF (TYPE.EQ.1) THEN
!!        -- gray scale
!         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
!      ELSE IF (TYPE.EQ.2) THEN
!!        -- rainbow
!         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
!      ELSE IF (TYPE.EQ.3) THEN
!!        -- heat
!         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
!      ELSE IF (TYPE.EQ.4) THEN
!!        -- weird IRAF
!         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
!      ELSE IF (TYPE.EQ.5) THEN
!!        -- AIPS
!         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
!      END IF
      ENDSUBROUTINE

endmodule galaxy
