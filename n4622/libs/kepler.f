!$Id: kepler.f 198 2011-09-21 02:45:04Z ccfeng $
        MODULE KEPLER
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8,SAVE                     ::eccen,angu_vel,rinit
        REAL*8,SAVE                     ::h = 10d-7
        REAL*8,SAVE                     ::G = 4.397d-6
        REAL*8,SAVE                     ::M = 1.08d11
        CONTAINS


        FUNCTION find_th(t)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                          ::k1,k2,t,th
        tt = 0.d0
        th = 0.d0
        Do while (tt .lt. t)
                k1 = h*iterating(th)
                k2 = h*iterating(th+k1)
                th = th + (k1 + k2)/2.d0
                tt = tt + h
                if(isNaN(k1)) then
                        write(*,*)'*********** th is nan '
                        write(*,*)th
                        write(*,*)'k1',k1,'k2',k2
                        STOP
                        exit
                ENDIF
        ENDDO
c       write(*,*)tt,th,radius(th)

        find_th = th
        END FUNCTION find_th

        FUNCTION iterating(th)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                  ::iterating,th
                iterating = (1.d0 + eccen*dcos(th))**2.d0
     c          /(1.d0 + eccen)**2.d0
                iterating = iterating*angu_vel
c               write(6,*)iterating
        END FUNCTION iterating

        FUNCTION radius(th)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                  ::radius,th
                radius = (1.d0 + eccen )/(1.d0 + eccen* dcos(th))*rinit
        END FUNCTION
        
        SUBROUTINE init_kepler(put_rinit,put_angular_velocity)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                  ::f
        INTEGER                 ::option = 1
        REAL*8                  ::confine

        !choose scenario
        if(option .eq. 0) then
                ! circle, decide by angular velocity
                angu_vel = put_angular_velocity
                rinit = (G*M*angu_vel**-2.d0)**(1.d0/3.d0)
        elseif(option .eq. 1) then
                ! circle, decide by initial position
                rinit = put_rinit
                !!!!!!decide retrograde or prograde HERE
                !!!!!!positive means prograde
                angu_vel = 1.d0*dsqrt(G*M*rinit**-3.d0)
        elseif(option .eq. 3) then
                ! decide all and check
                angu_vel = 1.d0*put_angular_velocity
                rinit = put_rinit
                confine = rinit**3.d0*angu_vel**2.d0
                if (confine .gt. 2.d0*G*M) then
                        write(6,*)'confine too large, stop'
                        write(6,*)'by angular, rinit max
     c                  is',(2.d0*G*M*angu_vel**-2.d0)**(1.d0/3.d0)
                        STOP
c               elseif(confine .le. G*M) then
c                       write(6,*)'confine too small, stop'
c                       write(6,*)'by angular, rinit min
c    c                  is',(G*M*angu_vel**-2.d0)**(1.d0/3.d0)

c                       STOP
                endif

                
        endif


        open(unit=70,file='./analysis/orbit_parameter')
        write(70,*)'angular velocity = ',angu_vel
        write(70,*)'rinit = ',rinit
        write(*,*)'angular velocity = ',angu_vel
        write(*,*)'rinit = ',rinit
c       f = find_focal()
c       eccen = f*(2.d0 - f)

        !!
        !!here
        !!
        f = (G*M - rinit**3.d0 * angu_vel**2.d0)
     c    / (rinit**2.d0*angu_vel**2.d0)

        eccen = f/(rinit+f)

        write(70,*)'eccentricity = ',eccen
        write(70,*)'focal length = ',f
        close(70)
        !!!!!!check eccentricity in range
        if((eccen > 1.0d0).or.(eccen < -10.d-15)) then
                write(6,*)'eccentricity out range'
                write(6,*)'eccentricity = ',eccen
                write(6,*)'focal = ',f
                stop
        endif

1000    open(unit=90,file='./analysis/orbit.dat')
        open(unit=80,file='./analysis/orbitxy.dat')
        

        END SUBROUTINE init_kepler

        SUBROUTINE print_kepler(t)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                          ::r,th
                r = radius(find_th(t))
                th = find_th(t)
                write(90,*)t,th,r
                write(80,*)t,r*dcos(th),r*dsin(th)
        END SUBROUTINE print_kepler

        FUNCTION output_angular_speed()
        IMPLICIT NONE
        REAL*8                          ::output_angular_speed
        output_angular_speed = angu_vel
        ENDFUNCTION output_angular_speed

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !               find focal length
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION find_focal()
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                          ::x

        x = 0.5d0
        write(6,*)'find focal start',focal(x,0)
        Do while(dabs(focal(x,0)) .ge. 10d-9) 
                x = x - focal(x,0)/focal(x,1)
                write(6,*)'find f',x
        enddo

        find_focal = x
        ENDFUNCTION find_focal

        FUNCTION focal(f,i)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        INTEGER                    ::i

        if (i .eq. 0) then
                focal = 
     -  f*G*M - 2*f**2*G*M + f**3*G*M + G*M*rinit - 
     -  2*f*G*M*rinit + f**2*G*M*rinit - 
     -  rinit**4*angu_vel**2
        elseif(i.eq.1)then
                 focal = 
     -  G*M - 4*f*G*M + 3*f**2*G*M - 2*G*M*rinit + 2*f*G*M*rinit
        elseif(i.eq.2)then
                 focal = 
     -  -4*G*M + 6*f*G*M + 2*G*M*rinit
        else
                write(*,*)'error'
        endif
        ENDFUNCTION focal

        SUBROUTINE end_kepler
        close(80)
        close(90)
        close(70)
        ENDSUBROUTINE

        END MODULE KEPLER

        MODULE companion_force
        CONTAINS

        FUNCTION companion_direct_force(x,y,mass_ratio,ieLL)
        USE kepler
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        REAL*8                  ::companion_direct_force(2)
        REAL*8                  ::mass_ratio

        gm    = (-1.d0)*g*m*mass_ratio
        th    = find_th(t)
        rad   = radius(th)
        cx    = rad*dcos(th)
        cy    = rad*dsin(th)

        roche = 1.26d0*10.d0/dble(REAL(ieLL))
        rsq   = x**2+y**2
          r   = dsqrt(rsq)
        cdx   = cx - x
        cdy   = cy - y
        drsq  = cdx**2.d0 + cdy**2.d0 + roche

        copf  = gm/drsq
        copfx = copf/dsqrt(drsq)*cdx
        copfy = copf/dsqrt(drsq)*cdy

        companion_direct_force = -(/copfx,copfy/)

        ENDFUNCTION companion_direct_force

        ENDMODULE companion_force
