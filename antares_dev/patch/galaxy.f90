        module galaxy
        DOUBLE PRECISION        ::SoundSpeed = 10.d0
        DOUBLE PRECISION        ::GravitatioalConst = 4.397d-6
        DOUBLE PRECISION        ::pi = 4.d0*datan(1.d0) 
        
        contains

        FUNCTION RC(r)
        IMPLICIT NONE
        DOUBLE PRECISION        ::RC,r

        RC = Vdisk(r) + Vhalo(r)

        ENDFUNCTION

        FUNCTION Vdisk(r)
        IMPLICIT NONE
        DOUBLE PRECISION        ::Vdisk,r
        DOUBLE PRECISION        ::G,M,c

        !Rotation Curve from Kuzmin disk
        !
        !    phi = -G M *(r^2+c^2)^-0.5
        !
        !    v   = sqrt(r*D[phi]) 
        !        
        !               G M R^2
        !        = (----------------)^0.5
        !            (c^2+r^2)^1.5
        !---------------------------

        G = GravitatioalConst
        M = 3.5d10
        c = 3.1d0

        Vdisk = dsqrt(R*(G*M/(R**2+c**2)**1.5))

        ENDFUNCTION

        FUNCTION Vhalo(r)
        IMPLICIT NONE
        DOUBLE PRECISION        ::Vhalo,r
        DOUBLE PRECISION        ::G,rho,a

        !Rotation Curve from Henquist halo
        !
        !           -4pi*G*rho*a^2
        !    phi = ------------------
        !              2(1+R/a)
        !
        !    v   = sqrt(r*D[phi]) 
        !        
        !              2pi G a^3 rho R
        !        = (---------------------)^0.5
        !                   (a+R)^2
        !---------------------------

        G   = GravitatioalConst
        rho = 4.11d7
        a   = 5.194d0

        Vhalo = dsqrt((2.d0*pi*G*a**3*rho*r)/(a+r)**2)

        ENDFUNCTION

        FUNCTION Linear_RC(r)
        IMPLICIT NONE
        DOUBLE PRECISION        ::Linear_RC,r

        DOUBLE PRECISION        ::v0 = 205.403
        DOUBLE PRECISION        ::e  = 0.0656488

        !Nearly flat Rotation Curve
        !
        !               r
        !    v = v0*(-------)^0.5
        !             r + e 
        !---------------------------

        Linear_RC =  v0 *(r/(r+e))**0.5

        ENDFUNCTION 

        FUNCTION density(r,z)
        IMPLICIT NONE
        DOUBLE PRECISION                ::density,r
        DOUBLE PRECISION,OPTIONAL       ::z

        density = GaussianDensity(r,z)
        !density = UniDensity(r,z)

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
        !  density = rho * exp(-r^2/(2a))
        !
        !---------------------------

        rho = 1.d0
!       a   = 7.d0
        a   = 10.d0

        GaussianDensity = rho * dexp(-r**2/2.d0/a)

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
