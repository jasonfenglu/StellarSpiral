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

rho = 19.33d0
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

endmodule galaxy
