        PROGRAM plot_rotation_curve
        USE rotation, only: angular,skapa
        IMPLICIT NONE
        REAL*8,ALLOCATABLE              ::x(:)
        REAL*8,ALLOCATABLE              ::y(:,:)
        INTEGER                         ::N = 1000
        INTEGER                         ::I
        REAL*8                          ::rang(4)
        REAL*8                          ::step
        REAL*8                          ::pspd = 4.19380495493055d0
        REAL*8,EXTERNAL                 ::find_res

        INTEGER                         ::pgopen
        CHARACTER(len=32)               ::arg
        CHARACTER(len=30)                    ::text

        !read argumet from input
        CALL get_command_argument(1, arg)
        IF (LEN_TRIM(arg) .ne. 0) then
                read(arg,*)pspd
        endif
                

        ALLOCATE(x(N))
        ALLOCATE(y(N,3))

        rang = (/0.,20.,0.,5./)
        step = (rang(2)-rang(1))/N
        DO I = 1, N
                x(I) = step*(I-1) + rang(1)
                y(I,1) = angular(x(I))
                y(I,2) = angular(x(I))-skapa(x(I),2)/2.d0
                y(I,3) = angular(x(I))+skapa(x(I),2)/2.d0
        ENDDO

        !setting draw env.
        if(PGOPEN('/XSERVE') .le. 0)stop
        CALL PGENV(0.1,20.,0.,100.,0,0)
        !plot rotation curve
        CALL PGLINE(N,REAL(x),REAL(y(:,1)))
        CALL PGSLS (2)
        CALL PGLINE(N,REAL(x),REAL(y(:,2)))
        CALL PGLINE(N,REAL(x),REAL(y(:,3)))

        CALL PGSLS(3)
        CALL PGSCI(5)
        y(:,3) = pspd
        CALL PGLINE(N,REAL(x),REAL(y(:,3)))
        write(text,'(a,g14.6)')'pattern speed is',pspd
        CALL PGTEXT(12.,80.,text)
        write(text,'(a,g14.6)')'resonance is at',find_res(pspd)
        CALL PGTEXT(12.,75.,text)

        DEALLOCATE(x)
        DEALLOCATE(y)
        write(*,*)'res is at', find_res(pspd),'for pattern speed:',pspd

        CALL PGCLOS
        STOP
        ENDPROGRAM plot_rotation_curve

        FUNCTION find_res(omega)
        IMPLICIT NONE
        REAL*8                  ::find_res,omega
        REAL*8                  ::B =  0.d0     !left bound
        REAL*8                  ::C = 20.d0     !right bound
        REAL*8                  ::R = 10.d0     !guess point
        REAL*8                  ::RE= 0.d0      !Relative error
        REAL*8                  ::AE= 0.d0      !Absolute error
        INTEGER                 ::IFLAG         !status

        CALL DFZERO(F,B,C,R,RE,AE,IFLAG)
        find_res = B

        CONTAINS 

        FUNCTION f(r)
        USE rotation, only: angular,skapa
        IMPLICIT NONE
        REAL*8                  ::f,r
        f = angular(r)-skapa(r,2)/2.d0-omega
        ENDFUNCTION


        ENDFUNCTION find_res
