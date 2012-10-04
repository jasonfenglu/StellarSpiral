        PROGRAM plot_smooth
        USE smooth
        IMPLICIT NONE
        REAL*8,ALLOCATABLE              ::x(:)
        REAL*8,ALLOCATABLE              ::y(:)
        INTEGER                         ::N = 1000
        INTEGER                         ::I
        REAL*8                          ::rang(4)
        REAL*8                          ::step

        INTEGER                         ::pgopen

        ALLOCATE(x(N))
        ALLOCATE(y(N))

        rang = (/0.,20.,0.,5./)
        step = (rang(2)-rang(1))/N
        DO I = 1, N
                x(I) = step*(I-1) + rang(1)
                y(I) = linear_increase(x(I),50.d0)
        ENDDO

        !setting draw env.
        if(PGOPEN('/xwin') .le. 0)stop
        CALL PGENV(0.1,20.,0.,1.2,0,0)
        CALL PGLINE(N,REAL(x),REAL(y))
        CALL PGCLOS

        DEALLOCATE(x)
        DEALLOCATE(y)

        STOP

        ENDPROGRAM plot_smooth
