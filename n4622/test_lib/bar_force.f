!$Id: companion_force.f 152 2011-08-04 07:49:59Z ccfeng $
        PROGRAM companion_potential_fourier
        USE force_from_companion_compoents
        USE kepler
        USE MATPLOT
        USE companion_force
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                    ::mass_ratio
        REAL*8,ALLOCATABLE              ::phi(:,:,:)
        INTEGER                         ::N1 = 30
        INTEGER                         ::N2 = 1000
        REAL*8                          ::rang(4)
        REAL*8                          ::force(2)
        INTEGER                         ::I,J
        REAL*8                          ::r

        INTEGER                         ::pgopen

        CHARACTER(3)                    ::mode



        rang = (/-20.,20.,-20.,20./)
        step = (rang(2)-rang(1))/REAL(N1)
        ALLOCATE(phi(N1,N1,2))

        Am   = 30.d0
        mass_ratio = 10.d-2
        t = 0.d0
        pspd = 20.d0
        !pragma omp parallel default(none)
        DO I = 1, N1
        !$OMP PARALLEL DO 
        DO J = 1, N1
                x = step*REAL(I-1) + rang(1)
                y = step*REAL(J-1) + rang(3) 

                 p02 = 1.d0

                  a2 = 5.2d0
              rsq    = x**2+y**2
                r    = dsqrt(rsq)
                p1   = rsq/(a2**2+rsq)**2
               psi   = p02*p1
               dpsi  = p02*2.d0*r*(a2**2-rsq)/(a2**2+rsq)**3
               cs2s  = (x**2-y**2)/rsq
               sn2s  = 2.d0*y*x/rsq
               cs2t  = cos(2.d0*pspd*t)
               sn2t  = sin(2.d0*pspd*t)
               cs2st = (cs2s*cs2t+sn2s*sn2t)
               sn2st = (sn2s*cs2t-cs2s*sn2t)
               px    = (dpsi*cs2st)*x/r+2.d0*psi*sn2st*y/rsq
               py    = (dpsi*cs2st)*y/r-2.d0*psi*sn2st*x/rsq
                
                force = (/px,py/)
                phi(I,J,:) = force
        ENDDO
        ENDDO
        call drawvec(REAL(phi(:,:,1)),REAL(phi(:,:,2)),N1,N1)
        DEALLOCATE(phi)
        STOP
        ENDPROGRAM
