      MODULE STELLARDISK
      DOUBLE PRECISION,PARAMETER::GravConst = 4.3d-6 
      DOUBLE PRECISION,PARAMETER::g = 4.3d0
      DOUBLE PRECISION,PARAMETER::pi=4.d0*atan(1.d0)
      DOUBLE PRECISION,SAVE     ::wr,wi
      type                      ::QParameterType
        DOUBLE PRECISION        ::Qod,q,rq
      ENDTYPE
      type(QParameterType),SAVE ::QParameter
      CONTAINS
      function Q(r)
      DOUBLE PRECISION  Q,r
      Qod = 1.d0
      q   = 6.2d0
      rq  = 2.5d0


      !Qod = QParameter%Qod
      !q   = QParameter%q  
      !rq  = QParameter%rq 

      Q = Qod*(1.d0 + q*dexp(-r**2/rq**2))
      endfunction

      function nu(r)
      IMPLICIT NONE
      DOUBLE COMPLEX    nu   
      DOUBLE PRECISION  r
      DOUBLE PRECISION  m 

      m = 2.d0

      nu = (dcmplx(wr,wi)-m*Omega(r))/kappa(r)
      
      endfunction

      function k3sqrt(r)
      IMPLICIT NONE
      DOUBLE COMPLEX            k3sqrt
      DOUBLE PRECISION,INTENT(in)::r
      DOUBLE PRECISION          ::rr
        
      k3sqrt = (dcmplx(kappa(r)/snsd(r)))**2*(dcmplx(Q(r))**-2 
     c         - 1.d0 + nu(r)**2)
!     rr = r +0.0000000001d0
!     k3sqrt = KappaOverASqr(rr)*(Q(rr)**(-2) - 1.d0 +
!    c nu(rr,wr,wi)**2+0.25d0*curF(rr)**2*q(rr)**2)

      endfunction
      
      function snsd(r)
      IMPLICIT NONE
      DOUBLE PRECISION          ::r,snsd
!     DOUBLE PRECISION,EXTERNAL ::Q


      snsd = Q(r)*pi*g*sigma0(r)/kappa(r)

      ENDFUNCTION

      function Sigma0(r)
      IMPLICIT NONE
      DOUBLE PRECISION          ::Sigma0,r
      DOUBLE PRECISION          ::m,a,b
      DOUBLE PRECISION          ::BOUND,EPSREL,EPSABS
      DOUBLE PRECISION          ::ans
      DOUBLE PRECISION          ::ABSERR
      INTEGER                   ::NEVAL,IERR,LIMIT,LENW,LAST,INF
      DOUBLE PRECISION,ALLOCATABLE ::WORK(:)
      INTEGER,ALLOCATABLE       ::IWORK(:)

      M     = 5.0d10
      a     = 2.7
      b     = 0.3

      BOUND = 0.d0
      INF   = 2
      EPSREL = 10d-10
      EPSABS = 10d-10
      LIMIT  = 100
      LENW   = LIMIT*4+2
      ALLOCATE(IWORK(LENW))
      ALLOCATE(WORK(LENW))
      CALL DQAGI(FUN,BOUND,INF,EPSABS,EPSREL,ANS,ABSERR,NEVAL,IERR,  
     c           LIMIT,LENW,LAST,IWORK,WORK)
      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)
      Sigma0 = ans/10d5
      contains 
      function FUN(z)
      IMPLICIT NONE
      DOUBLE PRECISION          ::fun,z
      fun =  
     c(b**2*M/4.d0/pi)* 
     c(a*r**2+(a+3.d0*sqrt(z**2+b**2))*(a+sqrt(z**2+b**2))**2)/ 
     c(r**2+(a+sqrt(z**2+b**2))**2)**2.5/(z**2+b**2)**1.5

      ENDFUNCTION

      ENDFUNCTION

      function kappa(r)
      IMPLICIT NONE
      DOUBLE PRECISION  kappa,r
      DOUBLE PRECISION  dr
      DOUBLE PRECISION  dOmega
!     DOUBLE PRECISION,EXTERNAL ::Omega
      dr = 0.00000000001d0
      dOmega = 0.d0
      dOmega = dOmega +  -3.d0/2.d0*Omega(r)
      dOmega = dOmega +        2.d0*Omega(r+dr)
      dOmega = dOmega +  -1.d0/2.d0*Omega(r+2*dr)

      kappa = sqrt(4.d0*Omega(r)**2*(1.d0+r/(2.d0*Omega(r))*dOmega))
      ENDFUNCTION

        function find_b(wr)
        IMPLICIT NONE
        DOUBLE PRECISION        find_b,wr
        DOUBLE PRECISION        ::B,C,i,RE,AE
!       DOUBLE PRECISION,EXTERNAL::kappa,Omega
        INTEGER                 ::IFLAG

        B = 1.d0
        C = 8.d0
        I = B
        

        CALL DFZERO(F,B,C,I,RE,AE,IFLAG)
        find_b = b

        CONTAINS

        FUNCTION F(r)
        DOUBLE PRECISION        ::F,r
        F = (2.d0*wr -2.d0*Omega(r))/kappa(r)

        F = F/kappa(r) -0.5d0
        ENDFUNCTION

        endfunction

        FUNCTION Omega(r)
        IMPLICIT NONE
        DOUBLE PRECISION          ::Omega,r
        !Halo
        DOUBLE PRECISION          ::Lh,rhoh,gHalo,VHalo
        !bulge
        DOUBLE PRECISION          ::rb,Mb,gBulge,VBulge
        !disk
        DOUBLE PRECISION          ::dM,da,db,VDisk
        !Halo
        Lh   = 0.3d0
        rhoh = 3.4e9
        gHalo = 4.d0*Lh**2.d0*pi*rhoh*(r - Lh*atan(r/Lh))
        gHalo = GravConst/(r**2)*gHalo
        VHalo = sqrt(r*gHalo)

        !Bulge
        Mb   = 12.0d8
        rb   = 0.8d0
        gBulge = 4*pi*rb**3*Mb*(-r/sqrt(1+r**2/rb**2)/rb+asinh(r/rb))
        gBulge = gBulge*GravConst/r**2
        VBulge = sqrt(r*gBulge)

        !Disk
        dM     = 5.0d10
        da     = 2.7
        db     = 0.3
        VDisk  = sqrt(dfunc(pDisk,r)*r)

        Omega  = sqrt(VHalo**2+VBulge**2+VDisk**2)/r
        CONTAINS
        FUNCTION pDisk(r)
        IMPLICIT NONE
        DOUBLE PRECISION        ::pDisk,r
        pDisk  = -GravConst*dM
        pDisk  = pDisk/sqrt(r**2+(da+db)**2)
        ENDFUNCTION

        ENDFUNCTION

        function dfunc(func,r)
        !
        ! Forward differential
        !
        IMPLICIT NONE
        DOUBLE PRECISION,EXTERNAL       ::func
        DOUBLE PRECISION                ::r,dfunc
        DOUBLE PRECISION                ::dr = 0.00000001d0

        dfunc = 0.d0
        dfunc = dfunc +  -3.d0/2.d0*func(r)
        dfunc = dfunc +        2.d0*func(r+dr)
        dfunc = dfunc +  -1.d0/2.d0*func(r+2*dr)
        dfunc = dfunc/dr
        endfunction

        ENDMODULE



!       function error(r,wr,wi)
!       IMPLICIT NONE
!       DOUBLE COMPLEX          ::error
!       DOUBLE COMPLEX,EXTERNAL ::k3sqrt
!       DOUBLE PRECISION,INTENT(in)::r,wr,wi
!       DOUBLE PRECISION        ::h=10d-5

!       error  = -(0.d0,1.d0)*sqrt(k3sqrt(r,wr,wi))
!       error  = error -
!    c  0.5d0/sqrt(k3sqrt(r,wr,wi))
!    c *(sqrt(k3sqrt(r+h,wr,wi))-sqrt(k3sqrt(r-h,wr,wi)))/(2.d0*h)

!       endfunction

!       FUNCTION curF(r)
!       IMPLICIT NONE
!       DOUBLE  PRECISION       ::curF,r
!       curF = 

!       ENDFUNCTION
