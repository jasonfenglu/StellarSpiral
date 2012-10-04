!$Id: companion_force.f 200 2011-09-21 09:09:58Z ccfeng $
        PROGRAM companion_potential_fourier
        USE force_from_companion_compoents
        USE kepler
        USE MATPLOT
        USE companion_force
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                    ::mass_ratio
        REAL*8,ALLOCATABLE              ::phi(:,:,:)
        INTEGER                         ::N1 = 30
        INTEGER                         ::N2 = 100
        REAL*8                          ::rang(4)
        REAL*8                          ::force(2)
        INTEGER                         ::I,J
        REAL*8                          ::r

        INTEGER                         ::pgopen

        CHARACTER(3)                    ::mode


        mode = '2'

        rang = (/-20.,20.,-20.,20./)
        step = (rang(2)-rang(1))/REAL(N1)
        ALLOCATE(phi(N1,N1,2))

        CALL init_kepler(30.d0,-2.d0)
        Am   = 30.d0
        mass_ratio = 10.d-2
        DO I = 1, N1
        DO J = 1, N1
                x = step*REAL(I-1) + rang(1)
                y = step*REAL(J-1) + rang(3) 
             force = force_from_companion(x,y,0.d0,Am,mass_ratio,mode)
             !force = companion_direct_force(x,y,mass_ratio,1024)
                phi(I,J,:) = force
        ENDDO
        ENDDO
        call drawvec(REAL(phi(:,:,1)),REAL(phi(:,:,2)),N1,N1)
        DEALLOCATE(phi)
        ALLOCATE(phi(N2,N2,2))
        step = (rang(2)-rang(1))/REAL(N2)
        DO I = 1, N2
        !$OMP DO PRIVATE(x,y,force,J,I)
        DO J = 1, N2
                x = step*REAL(I-1) + rang(1)
                y = step*REAL(J-1) + rang(3) 
             force = force_from_companion(x,y,0.d0,Am,mass_ratio,mode)
                phi(I,J,:) = force(:)
        ENDDO
        !$OMP END DO
        ENDDO

        DEALLOCATE(phi)
        STOP
        ENDPROGRAM
