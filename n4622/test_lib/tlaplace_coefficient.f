!$Id: tlaplace_coefficient.f 164 2011-08-15 07:04:32Z ccfeng $
        PROGRAM Laplace_coefficent
        USE slatec
        USE fourier_decom,fdlaplace_coefi=>dlaplace_coefi
     c  ,flaplace_coefi=>laplace_coefi
        IMPLICIT NONE
        REAL*8,ALLOCATABLE              ::x(:)
        REAL*8,ALLOCATABLE              ::y(:),yy(:)
        INTEGER                         ::N = 100
        INTEGER                         ::I,j
        REAL*8                          ::rang(4)
        REAL*8                          ::step

        INTEGER                         ::pgopen

        ALLOCATE(x(N))
        ALLOCATE(y(N))
        ALLOCATE(yy(N))

        rang = (/0.,1.,0.,5./)
        step = (rang(2)-rang(1))/N
        DO I = 1, N
                x(I) = step*(I-1) + rang(1)
                y(I) = dlaplace_coefi(2.d0,x(I))
                !y(I) = fdlaplace_coefi(2.d0,x(I))
                !y(I) = laplace_coefi1(x(I))
        ENDDO

        !setting draw env.
        !if(PGOPEN('/xwin') .le. 0)stop
        j = PGOPEN('/XSERVE')
        CALL PGENV(0.1,1.,-1.,7.,0,0)
        CALL PGLINE(N,REAL(x),REAL(y))

        DO I = 1, N
                x(I) = step*(I-1) + rang(1)
                yy(I) = fdlaplace_coefi(2.d0,x(I)) 
        ENDDO
        CALL PGLINE(N,REAL(x),REAL(yy))
        CALL PGCLOS
        write(*,*)y-yy
        DEALLOCATE(x)
        DEALLOCATE(y)
        DEALLOCATE(yy)

        STOP
        ENDPROGRAM
