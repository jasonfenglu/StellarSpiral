        module plotting
        CONTAINS
        SUBROUTINE plot2d(F,m,n,domain)
        IMPLICIT NONE
        DOUBLE PRECISION,ALLOCATABLE,INTENT(IN) ::F(:,:,:)!plotting data
        DOUBLE PRECISION                        ::domain(4)!plot range
        REAL                                    ::TR(6) !plot geometry
        REAL                                    ::vmax,vmin
        REAL                                    ::BRIGHT,CONTRA
        INTEGER                                 ::m,n   !dimentsion
        INTEGER                                 ::PGBEG


        IF (PGBEG(0,'/png',1,1) .NE. 1) STOP
        CALL PGSVP(0.0,0.95,0.0,0.95)


        TR(1) = REAL(dble(n)*domain(1)-domain(2))/REAL(n-1)
        TR(2) = REAL(domain(2)-domain(1))/REAL(n-1)
        TR(3) = 0.
        TR(4) = REAL(dble(m)*domain(3)-domain(4))/REAL(m-1)
        TR(5) = 0.
        TR(6) = REAL(domain(4)-domain(3))/REAL(m-1)

!       TR = (/0.,1.,0.,0.,0.,1./)
        vmax = 3.
        vmin = 0.
        

        BRIGHT = 0.5
        CONTRA = 2.0
        CALL PALETT(3,CONTRA,Bright)

        CALL PGBBUF
        CALL PGENV(real(domain(1)),real(domain(2)),real(domain(3)),real(domain(4)),0,0)
        CALL PGIMAG(REAL(F(:,:,3)),m,n,1,m,1,n,vmax,vmin,TR)
        CALL PGWEDG('RI', 1.0, 4.0, vmin, vmax, 'pixel value')
        CALL PGSCH(1.0)
        CALL PGEBUF

        CALL PGCLOS

        ENDSUBROUTINE

      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
!-----------------------------------------------------------------------
! Set a "palette" of colors in the range of color indices used by
! PGIMAG.
!-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
!
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
!
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
!
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
!
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
               0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
               0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
               0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
      IF (TYPE.EQ.1) THEN
!        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
!        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
!        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
!        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
!        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      ENDSUBROUTINE

        endmodule

PROGRAM find_all
USE PLOTTING
USE OMP_LIB
USE STELLARDISK
IMPLICIT NONE
type searchgrid_type
        DOUBLE PRECISION::coord(4,4,2)
        DOUBLE PRECISION::error(4,4)
endtype
type(searchgrid_type)            ::searchgrid
DOUBLE PRECISION                  ::wri,wii
INTEGER                          ::l

DOUBLE PRECISION,ALLOCATABLE    ::a(:,:)

!!$OMP PARALLEL DEFAULT(NONE),PRIVATE(a)
!!$OMP DO
!DO l = 1,4
!        ALLOCATE(a(2,2))
!        print *,'jer',l
!        DEALLOCATE(a)
!ENDDO
!!$OMP END DO
!!$OMP END PARALLEL
!stop


wri = 60.d0
wii = -0.5d0
CALL findpspsd(wri,wii)

STOP
!FUNCTION p(r)
!IMPLICIT NONE
!DOUBLE COMPLEX                  ::p
!DOUBLE PRECISION                ::r
!p = (0.d0,0.d0)
!ENDFUNCTION
!
!FUNCTION q(r)
!IMPLICIT NONE
!DOUBLE COMPLEX                  ::q
!DOUBLE COMPLEX,EXTERNAL         ::k3sqrt
!DOUBLE PRECISION                ::r
!q = k3sqrt(r,wr,wi)
!ENDFUNCTION
!
!SUBROUTINE search(wr,wi,err)
!IMPLICIT NONE
!DOUBLE COMPLEX,ALLOCATABLE      ::u(:,:),ui(:)
!DOUBLE COMPLEX,EXTERNAL         ::error
!DOUBLE PRECISION,EXTERNAL       ::find_b
!DOUBLE PRECISION                ::a,b
!DOUBLE PRECISION                ::h,r
!DOUBLE PRECISION,INTENT(IN)     ::wr,wi
!DOUBLE PRECISION                ::err
!INTEGER                         ::N=1000
!INTEGER                         ::i
!
!
!ALLOCATE(u(3,N))
!ALLOCATE(ui(3))
!a = 0.d0
!b = find_b(wr)
!ui = (/a,1.d0,0.d0/)
!h = (b-a)/REAL(n)
!CALL rk4(a,b,N,p,q,p,u,ui)
!err = dble(ABS(u(3,N)/u(2,N)-error(b,wr,wi)))
!write(*,*)wr,wi,err
!
!DEALLOCATE(u)
!DEALLOCATE(ui)
!
!
!ENDSUBROUTINE

END PROGRAM

SUBROUTINE single_grid(l,wri,wii)
USE STELLARDISK
IMPLICIT NONE
type searchgrid_type
        DOUBLE PRECISION::coord(12,12,2)
        DOUBLE PRECISION::error(12,12)
endtype
type(searchgrid_type)            ::searchgrid
DOUBLE PRECISION                ::dr,wri,wii,di
INTEGER                         ::l,i,j,p(2)


dr = 1.d0/10.0d0**(l-1)
di = 0.5d0/10.0d0**(l-1)
wri = wri +(-6.d0+0.5d0)*dr
wii = wii +(-6.d0+0.5d0)*dr
!CALL INIT_STELLARDISK(100,15.d0)
!print*,  abs(error())
!CALL ENDSTELLARDISK
!stop
DO i = 1,12
        searchgrid%coord(:,i,2) = dble(i-1)*dr + wii
        searchgrid%coord(i,:,1) = dble(i-1)*dr + wri
enddo

DO i = 1,12
DO j = 1,12
        wr = searchgrid%coord(i,j,1)
        wi = searchgrid%coord(i,j,2)
        CALL INIT_STELLARDISK(100,20.d0)
        searchgrid%error(i,j) = abs(error())
        CALL ENDSTELLARDISK
ENDDO
ENDDO

p = MINLOC(searchgrid%error(:,:))
i = p(1)
j = p(2)
wri = searchgrid%coord(i,j,1)
wii = searchgrid%coord(i,j,2)
!if(j.eq.1 .or. j.eq.12 .or. i.eq.1 .or. i.eq.12)l = l -1
!DO i = 1,12
!DO j = 1,12
!        print *,searchgrid%coord(i,j,:),searchgrid%error(i,j)
!ENDDO
!ENDDO
print *,wri,wii,searchgrid%error(i,j)

ENDSUBROUTINE

SUBROUTINE findpspsd(wri,wii)
IMPLICIT NONE
DOUBLE PRECISION                ::wri,wii
INTEGER                         ::l
do l = 1,5
        CALL single_grid(l,wri,wii)
enddo

ENDSUBROUTINE 
