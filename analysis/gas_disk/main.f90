        module plotting
        CONTAINS
        SUBROUTINE plot2d(F,m,n,x,y)
        IMPLICIT NONE
        DOUBLE PRECISION,ALLOCATABLE,INTENT(IN) ::F(:,:)!plotting data
        DOUBLE PRECISION                        ::x,y   !plot range
        REAL                                    ::TR(6) !plot geometry
        REAL                                    ::vmax,vmin
        INTEGER                                 ::m,n   !dimentsion
!       INTEGER                                 ::PGBEG


!     CALL PGMTXT('t',1.0,0.0,0.0,'Contouring using PGCONT')

      TR = (/0.,REAL(x/dble(m)),0.,0.,0.,REAL(y/dble(n))/)
!     TR = (/0.,1.,0.,0.,0.,1./)
        vmax = real(maxval(F))
        vmin = real(minval(F))
        
        write(*,*)'max,min',vmax,vmin

        CALL PGBBUF
        CALL PGENV(0.0,REAL(x),0.0,REAL(y),1,0)
        CALL PGGRAY(REAL(F),m,n,1,m,1,n,vmax,vmin,TR)
        CALL PGEBUF

        ENDSUBROUTINE

        endmodule

        module density
        DOUBLE PRECISION,SAVE                   ::r
        CONTAINS

        FUNCTION FUN(z)
        IMPLICIT NONE
        DOUBLE PRECISION                        ::FUN,z
        DOUBLE PRECISION                        ::G,M,a,b,Pi,Rho,h,c

        G = 4.3d-6
        M = 3.5d10
        a = 2.7d0
        b = 0.3d0
        rho=4.11d7
        h = 5.914
        c = 10.d0
        pi = 4.d0*datan(1.d0)
        
        FUN = &
        dexp(   &
        -1.d0/c**2* &
        (G*M/dsqrt((a+b)**2+r**2) - &
         G*M/dsqrt(r**2+(a+dsqrt(b**2+z**2))**2) + &
         2.d0*G*h**3*pi*rho/(h+r) - &
         2.d0*G*h**3*pi*rho/(h+dsqrt(r**2+z**2))  &
        ))

        ENDFUNCTION

        FUNCTION integral(rr)
        IMPLICIT NONE
        DOUBLE PRECISION                        ::integral
        DOUBLE PRECISION                        ::BOUND
        DOUBLE PRECISION                        ::EPSABS,EPSREL,ABSERR
        DOUBLE PRECISION                        ::ANS
        DOUBLE PRECISION,ALLOCATABLE            ::WORK(:)
        DOUBLE PRECISION                        ::rr
        INTEGER                                 ::NEVAL
        INTEGER                                 ::IERR,LIMIT,LENW,LAST
        INTEGER,ALLOCATABLE                     ::IWORK(:)
        INTEGER                                 ::INF


        
        r = rr
        BOUND  = 0.d0
        INF    = 2
        EPSABS = 1.d-15
        EPSREL = 1.d-15
        LIMIT  = 100
        LENW   = LIMIT*4+2
        ALLOCATE(WORK(LENW))
        ALLOCATE(IWORK(LENW))

        CALL DQAGI(FUN,BOUND,INF,EPSABS,EPSREL,ANS,ABSERR,NEVAL,IERR, &
                   LIMIT,LENW,LAST,IWORK,WORK)
!       write(*,*)ANS,ABSERR
        if(IERR.GT.0)write(*,*)'IERR',IERR,'!!!!!!!!!!!'
        if(IERR.eq.4)write(*,*)ABSERR
        integral = ANS
        DEALLOCATE(WORK)
        DEALLOCATE(IWORK)

        ENDFUNCTION

        endmodule

        MODULE FITTING
        type fit_model
                INTEGER                                 ::m!m eq.
                INTEGER                                 ::n!n unknows.
                DOUBLE PRECISION,ALLOCATABLE            ::datx(:)!data points
                DOUBLE PRECISION,ALLOCATABLE            ::daty(:) 
                DOUBLE PRECISION,ALLOCATABLE            ::x(:)!unknows
                DOUBLE PRECISION,ALLOCATABLE            ::fvec(:)!unknows
        endtype
        type(fit_model),save                            ::running_model

        CONTAINS
        SUBROUTINE FITTING_INIT(i,model,pmodel,m,datx,daty)
        IMPLICIT NONE
        type(fit_model)                         ::model
        type(fit_model)                         ::pmodel
        DOUBLE PRECISION                        ::datx(m)
        DOUBLE PRECISION                        ::daty(m)
        INTEGER                                 ::m
        INTEGER                                 ::i
        
        model%m = m
        model%n = 2

        ALLOCATE(model%datx(m))
        ALLOCATE(model%daty(m))
        model%datx = datx
        model%daty = daty

        ALLOCATE(model%x(model%n))
        if(i.gt.1)then
                model%x = pmodel%x
        else
                model%x = 2.d0
        endif
        ALLOCATE(model%fvec(m))

        running_model = model

        ENDSUBROUTINE

        SUBROUTINE FIT_CURVE(model)
        IMPLICIT NONE
        type(fit_model),INTENT(OUT)             ::model
        DOUBLE PRECISION                        ::tol
        INTEGER                                 ::info
        CALL lmdif1(fcn,running_model%m,running_model%n, &
                    running_model%x, running_model%fvec, tol,info)
        model = running_model
        ENDSUBROUTINE

        subroutine fcn ( m, n, x, fvec, iflag )
        IMPLICIT NONE
        DOUBLE PRECISION                        ::x(n)
        DOUBLE PRECISION                        ::fvec(m)
        INTEGER                                 ::iflag
        INTEGER                                 ::m,n

        fvec(:) = x(1) * dexp(-running_model%datx(:)**2/2.d0/x(2)**2) &
                - running_model%daty(:)

        ENDSUBROUTINE
        
        SUBROUTINE FIT_CLOS(model)
        IMPLICIT NONE
        type(fit_model)                         ::model
        DEALLOCATE(model%datx)
        DEALLOCATE(model%daty)
        DEALLOCATE(model%x)
        DEALLOCATE(model%fvec)

        ENDSUBROUTINE

        ENDMODULE

        PROGRAM disk
        use plotting
        use density,only:integral,FUN,fr=>r
        use fitting
        IMPLICIT NONE
        DOUBLE PRECISION,ALLOCATABLE            ::den(:,:)
        DOUBLE PRECISION,ALLOCATABLE            ::r(:),z(:)
        DOUBLE PRECISION                        ::lr,lz
        DOUBLE PRECISION                        ::rr,zz
        DOUBLE PRECISION                        ::dr,dz
        DOUBLE PRECISION                        ::integrate
        INTEGER                                 ::xnum,ynum
        INTEGER                                 ::i,j
        DOUBLE PRECISION                        ::b

        DOUBLE PRECISION                        ::xr(200),yr(200)
        INTEGER                                 ::PGBEG

        type(fit_model),ALLOCATABLE             ::model(:)

        xnum = 200
        ynum = 200

        lr = 20.d0
        lz = 3.d0

        dr = lr/dble(xnum)
        dz = lz/dble(ynum)

        b = 7.d0

        ALLOCATE(den(xnum,ynum))
        ALLOCATE(r(xnum))
        ALLOCATE(z(ynum))

        !find mid plane
        DO i = 1, xnum
                den(i,1) = integral(rr)**-1*dexp(-1.d0*rr**2/2.d0/b**2)
        ENDDo
        DO i = 1, xnum
        DO j = 2, ynum
                rr = dble(i) * dr
                zz = dble(j) * dz
                fr = rr
                den(i,j) = den(i,1)*FUN(zz)
        ENDDO
        ENDDO

        IF (PGBEG(0,'?',1,1) .NE. 1) STOP
        CALL PGSVP(0.0,0.95,0.0,0.95)

        !plot r-z density map
        CALL plot2d(den,xnum,ynum,lr,lz)

        !find scale length at different and plot
!       CALL PGENV(0.,REAL(lr),0.,REAL(lz),0,1)
        CALL PGBBUF
        ALLOCATE(model(xnum))
        DO J = 1,xnum
                xr = DBLE((/(I,I=1,ynum)/))*dz
                yr = den(J,:)
                CALL FITTING_INIT(J,model(J),model(J-1),ynum,xr,yr)
                CALL FIT_CURVE(model(J))
        ENDDO

        xr = DBLE((/(I,I=1,xnum)/))*dr
        DO J = 1,xnum
                yr(J) = model(J)%x(2)
        ENDDO
        CALL PGSCI(3)
        CALL PGPT(xnum,REAL(xr),REAL(yr),-3)
        CALL PGEBUF
        CALL PGSCI(1)
        CALL PGENV(0.,REAL(lr),0.,REAL(lz),0,1)
        CALL PGPT(xnum,REAL(xr),REAL(yr),-3)
        CALL PGCLOS
        print *,model(1)%x(2)*10.d2,model(1)%x(1)

1000    DEALLOCATE(den)
        DO J = 1,xnum
        CALL FIT_CLOS(model(J))
        ENDDO
        STOP
        END