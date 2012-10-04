!$Id: companion_potential_fourier.f 79 2011-07-12 04:38:26Z ccfeng $
        PROGRAM companion_potential_fourier
        USE potential_from_companion
        USE kepler
        USE MATPLOT
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                    ::mass_ratio
        REAL*8,ALLOCATABLE              ::phi(:,:,:)
        INTEGER                         ::N1 = 500
        INTEGER                         ::N2 = 1000
        REAL*8                          ::rang(4)
        REAL*8                          ::force(2)
        INTEGER                         ::I,J
        REAL*8                          ::r

        INTEGER                         ::pgopen


        rang = (/-20.,20.,-20.,20./)
        step = (rang(2)-rang(1))/REAL(N1)
        ALLOCATE(phi(N1,N1,2))

        CALL init_kepler(30.d0,-2.d0)
        Am   = 30.d0
        mass_ratio = 10.d-2
        !pragma omp parallel default(none)
        DO I = 1, N1
        !$OMP PARALLEL DO 
        DO J = 1, N1
                x = step*REAL(I-1) + rang(1)
                y = step*REAL(J-1) + rang(3) 
                phi(I,J,1) = potential_1(x,y,0.d0,Am)
                phi(I,J,2) = potential_i(x,y,0.d0,Am,mass_ratio)
        ENDDO
        ENDDO
        CALL drawmat(REAL(phi(:,:,1)),N1,N1)
        CALL drawmat(REAL(phi(:,:,2)),N1,N1)
        CALL drawmat(REAL(phi(:,:,1)+phi(:,:,2)),N1,N1)
        STOP
        ENDPROGRAM
