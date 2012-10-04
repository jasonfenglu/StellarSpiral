!$Id: rotation_data.f 111 2011-07-20 03:12:57Z ccfeng $
      MODULE rotation
      type rotparameters_bulge
              REAL*8                    ::v0 = 222.841d0
              REAL*8                    ::A  = -0.0205608d0
              REAL*8                    ::B  = 0.268375d0
      end type rotparameters_bulge

      type rotparameters_old  
              REAL*8                    ::v0 = 108.627d0
              REAL*8                    ::A  = 0.590664d0
              REAL*8                    ::B  = 0.406462d0
      end type rotparameters_old

      REAL*8,SAVE                       ::v0
      REAL*8,SAVE                       ::A 
      REAL*8,SAVE                       ::B 
      CONTAINS


ccccc Functions:ccccccccccccccccccccccccccccccccccccc
      FUNCTION rotcurve(r)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      !choose rotation curve here
      type(rotparameters_bulge) ::parameters
      !init rotation curve parameters
              v0 = parameters%v0
              A  = parameters%A
              B  = parameters%B
      !output rotcurve
      rotcurve = r**(A+1.d0)*v0/(r+r**(A+B))
      END FUNCTION rotcurve

      FUNCTION angular(r)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      angular = rotcurve(r)/r
      END FUNCTION angular

      FUNCTION skapa(r,m)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      t1 = r**A*dSqrt(2.d0*(1.d0 + A)*r - 2.d0*(-2.d0 + B)*r**(A + B))
     &*v0
      t2 = (r + r**(A + B))**1.5
      skapa = t1/t2
      
      END FUNCTION skapa

      FUNCTION pspd(r,m)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      select case(m)
      case (1)
             pspd = angular(r) - skapa(r,m)
      case (2)
            pspd = angular(r) + skapa(r,m)/m
      case default
            write(*,*)'pattern speed error!!'
      end select

      END FUNCTION pspd
      END MODULE rotation
