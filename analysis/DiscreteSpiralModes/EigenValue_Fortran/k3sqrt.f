      MODULE STELLARDISK
      CONTAINS
      function Q(r)
      DOUBLE PRECISION  Q,r
      Qod = 1.d0
      q   = 3.9d0
      rq  = 3.1d0
      Q = Qod*(1.d0 + q*dexp(-r**2/rq**2))
      endfunction

      function nu(r,wr,wi)
      IMPLICIT NONE
      DOUBLE COMPLEX    nu   
      DOUBLE PRECISION  wr,wi,r
      DOUBLE PRECISION  m 

      m = 2.d0

      nu = (dcmplx(wr,wi)-m*Omega(r))/kappa(r)
      
      endfunction

      function k3sqrt(r,wr,wi)
      IMPLICIT NONE
      DOUBLE COMPLEX            k3sqrt
!     DOUBLE COMPLEX,EXTERNAL   ::nu
      DOUBLE PRECISION,INTENT(in)::r,wr,wi
      DOUBLE PRECISION          ::rr
!     DOUBLE PRECISION,EXTERNAL ::KappaOverASqr,Q,curF
        
      k3sqrt = (kappa(r)/snsd(r))**2*(Q(r)**-2 - 1.d0 + nu(r,wr,wi)**2)
!     rr = r +0.0000000001d0
!     k3sqrt = KappaOverASqr(rr)*(Q(rr)**(-2) - 1.d0 +
!    c nu(rr,wr,wi)**2+0.25d0*curF(rr)**2*q(rr)**2)


      endfunction
      
      function snsd(r)
      IMPLICIT NONE
      DOUBLE PRECISION          ::r,snsd
!     DOUBLE PRECISION,EXTERNAL ::Q
      DOUBLE PRECISION          ::pi=4.d0*atan(1.d0)
      DOUBLE PRECISION          ::g = 4.3


      snsd = Q(r)*pi*g*sigma0(r)/kappa(r)

      ENDFUNCTION


      function Sigma0(r)
      IMPLICIT NONE
      DOUBLE PRECISION          ::Sigma0,r
      DOUBLE PRECISION          ::ssigma0,rd
      ssigma0 = 1206.d0
      rd      = 2.83d0

      Sigma0 = ssigma0*dexp(-r/rd)
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

        function error(r,wr,wi)
        IMPLICIT NONE
        DOUBLE COMPLEX          ::error
!       DOUBLE COMPLEX,EXTERNAL ::k3sqrt
        DOUBLE PRECISION,INTENT(in)::r,wr,wi
        DOUBLE PRECISION        ::h=10d-5

        error  = -(0.d0,1.d0)*sqrt(k3sqrt(r,wr,wi))
        error  = error -
     c  0.5d0/sqrt(k3sqrt(r,wr,wi))
     c *(sqrt(k3sqrt(r+h,wr,wi))-sqrt(k3sqrt(r-h,wr,wi)))/(2.d0*h)

        endfunction

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
        DOUBLE PRECISION          ::GravConst = 4.3d-6 
        DOUBLE PRECISION          ::Lh,rhoh,gHalo,VHalo
        DOUBLE PRECISION          ::Mb,rm,gBulge,VBulge
        DOUBLE PRECISION          ::ssigma0,Rd,sigma0,y,VDisk
        DOUBLE PRECISION          ::pi,g = 4.3d0
        DOUBLE PRECISION,EXTERNAL ::bessel_k,bessel_i
        pi = 4.d0*atan(1.d0)
        !Halo
        Lh   = 3.d0;
        rhoh = 3.3e7;
        gHalo = 4.d0*Lh**2.d0*pi*rhoh*(r - Lh*atan(r/Lh));
        gHalo = GravConst/(r**2)*gHalo;
        VHalo = sqrt(r*gHalo);

        !Bulge
        Mb   = 1.5e10;
        rm   = 1.5;
        gBulge = Mb*r/(r+rm);
        gBulge = GravConst/(r**2)*gBulge;
        VBulge = sqrt(r*gBulge);

        !Disk
        ssigma0 = 1206.d0;
        Rd     = 2.83;
        sigma0 = Sigma0*exp(-r/Rd);
        y      = r/2.d0/Rd;
        VDisk  = 4.d0*g*Sigma0*Rd*y**2.d0
        VDisk  = VDisk*(bessel_i(0,y)*bessel_k(0,y)-
     c           bessel_i(1,y)*bessel_k(1,y));
        VDisk  = sqrt(VDisk);

        Omega  = sqrt(VHalo**2+VBulge**2+VDisk**2)

        ENDFUNCTION

!       FUNCTION curF(r)
!       IMPLICIT NONE
!       DOUBLE  PRECISION       ::curF,r
!       curF = 

!       ENDFUNCTION
        ENDMODULE
