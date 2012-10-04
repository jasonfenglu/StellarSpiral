!$Id: companion_forcing_fourier_decompose.f 166 2011-08-15 10:12:27Z ccfeng $
        MODULE fourier_decom
        REAL*8,PARAMETER::pi = datan(1.d0)*4
        CONTAINS

        FUNCTION laplace_coefi0(beta)
        USE MATH
        IMPLICIT NONE
	REAL*8          ::beta
	REAL*8          ::laplace_coefi0

        laplace_coefi0 = 4.d0/pi/(1.d0+beta)
     &  *ellipticK(4.d0*beta/(1.d0+beta)**2.d0)

        ENDFUNCTION laplace_coefi0

        FUNCTION dlaplace_coefi0(beta)
        USE MATH
        IMPLICIT NONE
        REAL*8          ::beta
        REAL*8          ::dlaplace_coefi0
        REAL*8          ::phase

        phase = 4.d0*beta/(1.d0+beta)**2

        dlaplace_coefi0 = 
     &  ((1.d0+beta)*ellipticE(phase)+(-1.d0+beta)*ellipticK(phase))
     &  *(-2.d0)/(pi*beta*(beta**2-1.d0))

       
        ENDFUNCTION dlaplace_coefi0

        FUNCTION laplace_coefi1(beta)
        USE MATH
        IMPLICIT NONE
        REAL*8          ::beta
        REAL*8          ::laplace_coefi1
        REAL*8          ::phase

        phase = 4.d0*beta/(1.d0+beta)**2

        laplace_coefi1 = 
     & ((1.d0+beta**2)*ellipticK(phase)-(1.d0+beta)**2*ellipticE(phase))
     &  *2.d0/pi/beta/(1.d0+beta)
        ENDFUNCTION laplace_coefi1

        FUNCTION dlaplace_coefi1(beta)
        USE MATH
        IMPLICIT NONE
        REAL*8          ::beta
        REAL*8          ::dlaplace_coefi1
        REAL*8          ::phase

        phase = 4.d0*beta/(1.d0+beta)**2

        dlaplace_coefi1 = 
     &  ((1.d0+beta)*ellipticE(phase)+(-1.d0+beta)*ellipticK(phase))
     &  *(-2.d0)/(pi*beta**2*(beta**2-1.d0))

!       dlaplace_coefi1 = 
!    &  2.d0*(beta+1.d0)/pi/beta
!    &  *(ellipticK(phase) - ellipticE(phase))
!       
        ENDFUNCTION dlaplace_coefi1

        FUNCTION laplace_coefi(m,beta)
!DEC$ ATTRIBUTES ALIAS : 'wrapper_laplace_coei' :: wrapper_laplace_coei
        IMPLICIT NONE
        REAL*8                  ::laplace_coefi
        REAL*8                  ::m,beta,ans
                call wrapper_laplace_coei(m,beta,ans)
                laplace_coefi = ans
        ENDFUNCTION laplace_coefi

        FUNCTION dlaplace_coefi(m,beta)
!DEC$ ATTRIBUTES ALIAS:'wrapper_dlaplace_coei'::wrapper_dlaplace_coei
        IMPLICIT NONE
        REAL*8                  ::dlaplace_coefi
        REAL*8                  ::m,beta,ans
                call wrapper_dlaplace_coei(m,beta,ans)
                dlaplace_coefi = ans
        ENDFUNCTION dlaplace_coefi

        ENDMODULE fourier_decom

        MODULE force_from_companion_compoents
        CONTAINS

        FUNCTION force_from_companion(x,y,t,Am,mass_ratio,option)
        USE kepler
        USE fourier_decom,flaplace_coefi=>laplace_coefi
     c                  ,fdlaplace_coefi=>dlaplace_coefi
        USE slatec,only:laplace_coefi,dlaplace_coefi
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                    ::force_from_companion(2)
        REAL*8                    ::mass_ratio
        REAL*8                    ::runit(2)
        REAL*8                    ::thunit(2)
        REAL*8                    ::force(2,4) = 0.d0
        REAL*8                    ::beta
        CHARACTER(*)              ::option


        gm       = G*M*mass_ratio
        r        = dsqrt(x**2 + y**2)
        th       = datan2(y,x)
        runit    = (/x/r,y/r/)
        thunit   = (/-y/r,x/r/)
        beta     = r/Am
        bbeta    = 4.d0*beta/(1.d0+beta)**2

        if (scan(option,'0')/=0)then
                dphi0    = 5.d-1*dlaplace_coefi0(beta)*gm/(Am**2)
                force(:,1) = dphi0*runit
        endif

        if (scan(option,'1')/=0)then
        dphi1_r  = dlaplace_coefi1(beta)*gm/(Am**2)
     &  *dcos(angu_vel*t-th)
        dphi1_th = laplace_coefi1(beta)*gm/(Am*r)
     &  *dsin(angu_vel*t-th) 
        force(:,2) =  dphi1_r*runit + dphi1_th*thunit
        endif

        if (scan(option,'2')/=0)then
        dphi2_r  = dlaplace_coefi(2.d0,beta)*gm/(Am**2)
     &  *dcos(2.d0*(angu_vel*t-th))
        dphi2_th = laplace_coefi(2.d0,beta)*gm/(Am*r)
     &  *2.d0*dsin(2.d0*(angu_vel*t-th))
        force(:,3) =  dphi2_r*runit + dphi2_th*thunit
        endif

        if (scan(option,'I')/=0)then
        dphi_inde_r = -mass_ratio*Am*angu_vel**2
     &  *dcos(angu_vel*t-th)
        dphi_inde_th= -mass_ratio*Am*angu_vel**2
     &  *dsin(angu_vel*t-th)
        force(:,3) = dphi_inde_r*runit + dphi_inde_th*thunit
        endif

        force_from_companion = SUM(force,2)

        ENDFUNCTION force_from_companion

        ENDMODULE force_from_companion_compoents

        MODULE potential_from_companion
        CONTAINS

        FUNCTION potential_1(x,y,t,Am)
        USE fourier_decom
        USE kepler
        IMPLICIT NONE
        REAL*8                  ::potential_1
        REAL*8                  ::x,y,r,th,t
        REAL*8                  ::beta,Am

        r = dsqrt(x**2 + y**2)
        th= datan2(y,x)
        beta = r/Am

        potential_1 = laplace_coefi1(beta) 
     &  * dcos(angu_vel*t-th)*(-g*m/Am)
        ENDFUNCTION potential_1
        
        FUNCTION potential_i(x,y,t,Am,mass_ratio)
        USE fourier_decom
        USE kepler
        IMPLICIT NONE
        REAL*8                  ::potential_i
        REAL*8                  ::x,y,r,th,t
        REAL*8                  ::beta,Am,mass_ratio

        r = dsqrt(x**2 + y**2)
        th= datan2(y,x)
        beta = r/Am

        potential_i = mass_ratio*Am*angu_vel**2*r
     &  * dcos(angu_vel*t-th)
        ENDFUNCTION potential_i

        ENDMODULE potential_from_companion
