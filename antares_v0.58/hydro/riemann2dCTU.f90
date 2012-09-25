subroutine riemann2dCTU(q_loc)
use common_params
implicit none
integer::i,j,k
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision,dimension(:,:,:),allocatable::DQ,GC
double precision,dimension(:,:),allocatable::APDQ,AMDQ,FC
double precision,dimension(:,:)  ,allocatable::Q1D
double precision::dtddx,dtddy
double precision::ql(NVAR),qr(NVAR)
integer::NBC,NX1D


allocate(DQ(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR))

DQ=0.d0
dtddx = dt/dx
dtddy = dt/dy

!!! x-sweeps
NBC=ibuf
NX1D=ncell_loc(1)
allocate( Q1D(1-NBC:NX1D+NBC,NVAR))
allocate(APDQ(2-NBC:NX1D+NBC,NVAR))
allocate(AMDQ(2-NBC:NX1D+NBC,NVAR))
allocate(  FC(2-NBC:NX1D+NBC,NVAR))
allocate(  GC(1-NBC:NX1D+NBC,NVAR,2))

APDQ=0.d0
AMDQ=0.d0
FC=0.d0
GC=0.d0

do j=0,ncell_loc(2)+1
!!! copy information from each slice
  do i=1-ibuf,ncell_loc(1)+ibuf
    do k=1,NVAR
     Q1D(i,k)=q_loc(i,j,k)
    enddo
  enddo

  call flux2dCTU(NX1D,NBC,1,Q1D,AMDQ,APDQ,FC,GC)
  
  do i=1,ncell_loc(1)
    do k=1,NVAR
     DQ(i,j,k)=DQ(i,j,k) &
              -dtddx*((APDQ(i,k)+AMDQ(i+1,k))+ &
                      (FC(i+1,k)-FC(i,k))) &
              -dtddy*(GC(i,k,2)-GC(i,k,1))
     DQ(i,j-1,k)=DQ(i,j-1,k)-dtddy*GC(i,k,1)
     DQ(i,j+1,k)=DQ(i,j+1,k)+dtddy*GC(i,k,2)
    enddo 
  enddo
enddo
deallocate(Q1D)
deallocate(APDQ)
deallocate(AMDQ)
deallocate(FC)
deallocate(GC)

!!! y-sweeps

NBC=jbuf
NX1D=ncell_loc(2)
allocate( Q1D(1-NBC:NX1D+NBC,NVAR))
allocate(APDQ(2-NBC:NX1D+NBC,NVAR))
allocate(AMDQ(2-NBC:NX1D+NBC,NVAR))
allocate(  FC(2-NBC:NX1D+NBC,NVAR))
allocate(  GC(1-NBC:NX1D+NBC,NVAR,2))
APDQ=0.d0
AMDQ=0.d0
FC=0.d0
GC=0.d0

do i=0,ncell_loc(1)+1
!!! copy information from each slice
  do j=1-jbuf,ncell_loc(2)+jbuf
    do k=1,NVAR
     Q1D(j,k)=q_loc(i,j,k)
    enddo
  enddo

  call flux2dCTU(NX1D,NBC,2,Q1D,AMDQ,APDQ,FC,GC)

  do j=1,ncell_loc(2)
    do k=1,NVAR
     DQ(i,j,k)=DQ(i,j,k) &
              -dtddy*((APDQ(j,k)+AMDQ(j+1,k))+ &
                     (FC(j+1,k)-FC(j,k))) &
              -dtddx*(GC(j,k,2)-GC(j,k,1))
     DQ(i-1,j,k)=DQ(i-1,j,k)-dtddx*GC(j,k,1)
     DQ(i+1,j,k)=DQ(i+1,j,k)+dtddx*GC(j,k,2)
    enddo
  enddo
enddo
deallocate(Q1D)
deallocate(APDQ)
deallocate(AMDQ)
deallocate(FC)
deallocate(GC)
!!!!!!!!!!!!
do k=1,3
  do j=1,ncell_loc(2)
    do i=1,ncell_loc(1)
      q_loc(i,j,k)=q_loc(i,j,k)+DQ(i,j,k)
    enddo
  enddo     
enddo


deallocate(DQ)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine flux2dCTU(NX1D,NBC,IXY,Q1D,AMDQ,APDQ,FC,GC)
use common_params
implicit none
integer::IXY,NW,NX1D,I,L,K,NBC
double precision:: Q1D(1-NBC:NX1D+NBC,NVAR)
double precision::AMDQ(2-NBC:NX1D+NBC,NVAR)
double precision::APDQ(2-NBC:NX1D+NBC,NVAR)
double precision::FC(2-NBC:NX1D+NBC,NVAR)
double precision::GC(1-NBC:NX1D+NBC,NVAR,2)
double precision::TEMP,DTDX1D
double precision,dimension(:,:),allocatable::S
double precision,dimension(:,:,:),allocatable::W
double precision,dimension(:),allocatable::US,VS
double precision,allocatable,dimension(:,:)::BMASDQ
double precision,allocatable,dimension(:,:)::BPASDQ

NW=NVAR
allocate(S(2-NBC:NX1D+NBC,NW))
allocate(W(2-NBC:NX1D+NBC,NVAR,NW))
allocate(US(2-NBC:NX1D+NBC))
allocate(VS(2-NBC:NX1D+NBC))
allocate(BMASDQ(2-NBC:NX1D+NBC,NVAR))
allocate(BPASDQ(2-NBC:NX1D+NBC,NVAR))

if(IXY .eq. 1) then
  DTDX1D=dt/dx
elseif(IXY .eq. 2) then
  DTDX1D=dt/dy
endif

FC=0.d0
GC=0.d0

call RPN2(NX1D,NBC,NW,IXY,Q1D,AMDQ,APDQ,S,W,US,VS)
call LIMIT2(W,NX1D,NBC,NW,S)

do I=2-NBC,NX1D+NBC
  do L=1,NW
    TEMP=0.5d0*dabs(S(I,L))*(1.d0-DTDX1D*dabs(S(I,L)))
    do K=1,NW
      FC(I,K)=FC(I,K)+TEMP*W(I,K,L)
    enddo
  enddo
enddo

call RPT2(NX1D,NBC,IXY,1,Q1D,AMDQ,BMASDQ,BPASDQ,US,VS)

do I=2-NBC,NX1D+NBC
   do K=1,NW
     GC(I-1,K,1)=GC(I-1,K,1)-BMASDQ(I,K)
     GC(I-1,K,2)=GC(I-1,K,2)-BPASDQ(I,K)
   enddo
enddo

call RPT2(NX1D,NBC,IXY,2,Q1D,APDQ,BMASDQ,BPASDQ,US,VS)
do I=2-NBC,NX1D+NBC
   do K=1,NW
     GC(I,K,1)=GC(I,K,1)-BMASDQ(I,K)
     GC(I,K,2)=GC(i,K,2)-BPASDQ(I,K)
   enddo
enddo

TEMP=0.5d0*DTDX1D
do I=1-NBC,NX1D+NBC
   do K=1,NW
     GC(I,K,1)=TEMP*GC(I,K,1)
     GC(I,K,2)=TEMP*GC(I,K,2)
   enddo
enddo

deallocate(S)
deallocate(W)
deallocate(US)
deallocate(VS)
deallocate(BMASDQ)
deallocate(BPASDQ)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Modified From Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RPN2(NX1D,NBC,NW,IXY,Q1D,AMDQ,APDQ,S,W,US,VS)
use common_params
implicit none
integer::NX1D,NBC,IXY,NW
double precision::Q1D(1-NBC:NX1D+NBC,NVAR)
double precision::AMDQ(2-NBC:NX1D+NBC,NVAR)
double precision::APDQ(2-NBC:NX1D+NBC,NVAR)
double precision::S(2-NBC:NX1D+NBC,NW),W(2-NBC:NX1D+NBC,NVAR,NW)
double precision::US(2-NBC:NX1D+NBC),VS(2-NBC:NX1D+NBC)

double precision:: TOL,EXPMAX
double precision::F1,F2,F3,RHO,U,V
double precision::RHOR,UR,VR
double precision::RHOL,UL,VL
double precision::RHOM,UM,RHOS,EPS
double precision::C1,C2
double precision::PHIL,PHIR
integer::IU,IV
logical::ERR
integer::I,K

F1(RHO,U,V)=RHO*U
F2(RHO,U,V)=RHO*(U*U+snd*snd)
F3(RHO,U,V)=RHO*U*V


TOL = 1.d-12
EXPMAX=690.d0
! Check the normal direction.
      IF (IXY .EQ. 1) THEN
         IU = 2
         IV = 3
      ELSE IF (IXY .EQ. 2) THEN
         IU = 3
         IV = 2
      ELSE
         PAUSE 'Unknown sweep code.'
      END IF

!C Solve the Riemann problem at each interface.
!      SMAX = 0.0D0
      RHOR = Q1D(1-NBC, 1)
        UR = Q1D(1-NBC, IU) / RHOR
        VR = Q1D(1-NBC, IV) / RHOR

      DO I = 2 - NBC, NX1D + NBC
!C        Find the middle state.
         RHOL = RHOR
         UL = UR
         VL = VR
         RHOR = Q1D(I,1)
         UR = Q1D(I,IU) / RHOR
         VR = Q1D(I,IV) / RHOR
!C        All-shock solution
         C1 = DSQRT(RHOL * RHOR)
         C2 = (UL - UR) * C1 / (snd * (DSQRT(RHOL) + DSQRT(RHOR)))
         RHOM = 0.25D0 * (C2 + DSQRT(C2 * C2 + 4.0D0 * C1)) ** 2
         IF (RHOM .LT. RHOL .OR. RHOM .LT. RHOR) THEN
!C           All-rarefaction solution
            C1 = (UL - UR) / snd 
            IF (C1 .GT. EXPMAX) PAUSE 'Difficulty in all-rarefaction solution.'
            RHOM = DSQRT(RHOL * RHOR * DEXP(C1))
            IF (RHOM .GT. RHOL .OR. RHOM .GT. RHOR) THEN
!C              Otherwise
               CALL NEWTON(DMIN1(RHOL,RHOR),TOL,RHOM,ERR,RHOL,UL,VL,RHOR,UR,VR)
               IF (ERR) then
                 write(*,*) 'Newton''s method failed.'
                 stop
               endif
            END IF
         END IF
         UM = 0.5D0 * (PHIL(RHOM,RHOL,UL,VR) + PHIR(RHOM,RHOR,UR,VR))

!C        Wave speeds
         IF (RHOM .GT. RHOL) THEN
!C           1-shock
            EPS = (RHOM - RHOL) / RHOL
            IF (EPS .GT. TOL) THEN
               S(I,1) = (RHOM * UM - RHOL * UL) / (RHOM - RHOL)
            ELSE
!C              Weak shock
               S(I,1) = UL - snd * (1.0D0 + 0.5D0 * EPS &
                                       - 0.125D0 * EPS * EPS)
            END IF
         ELSE
!C           1-rarefaction
            S(I,1) = 0.5D0 * (UM + UL) - snd 
         END IF
         S(I,2) = UM
         IF (RHOM .GT. RHOR) THEN
!C           3-shock
            EPS = (RHOM - RHOR) / RHOR
            IF (EPS .GT. TOL) THEN
               S(I,3) = (RHOM * UM - RHOR * UR) / (RHOM - RHOR)
            ELSE
!C              Weak shock
               S(I,3) = UR + snd * (1.0D0 + 0.5D0 * EPS  &
                                       - 0.125D0 * EPS * EPS)
            END IF
         ELSE
!C           3-rarefaction
            S(I,3) = 0.5D0 * (UM + UR) + snd
         END IF
  !       DO 20 K = 1, NW
            !IF (DABS(S(I,K)) .GT. SMAX) SMAX = DABS(S(I,K))
  ! 20    CONTINUE
!C        Waves
         W(I,1, 1) = RHOM - RHOL
         W(I,IU,1) = RHOM * UM - RHOL * UL
         W(I,IV,1) = W(I,1,1) * VL
         W(I,1, 2) = 0.0D0
         W(I,IU,2) = 0.0D0
         W(I,IV,2) = RHOM * (VR - VL)
         W(I,1, 3) = RHOR - RHOM
         W(I,IU,3) = RHOR * UR - RHOM * UM
         W(I,IV,3) = W(I,1,3) * VR

!C        Find the state with zero speed.
         IF ((RHOM .GT. RHOL .AND. S(I,1) .GE. 0.0D0) .OR. &
            (RHOM .LE. RHOL .AND. UL - snd .GE. 0.0D0)) THEN
!C           Left state
            RHOS = RHOL
            US(I) = UL
            VS(I) = VL
         ELSE IF (RHOM .LE. RHOL .AND. UM - snd .GE. 0.0D0) THEN
!C           1-centered rarefaction wave
            RHOS = RHOL * DEXP((UL - snd) / snd)
            US(I) = snd 
            VS(I) = VL
         ELSE IF ((RHOM .GT. RHOR .AND. S(I,3) .GE. 0.0D0) .OR.  &
                 (RHOM .LE. RHOR .AND. UM + snd .GE. 0.0D0)) THEN
!C           Middle state
            RHOS = RHOM
            US(I) = UM
            IF (UM .GE. 0.0D0) THEN
               VS(I) = VL
            ELSE
               VS(I) = VR
            END IF
         ELSE IF (RHOM .LE. RHOR .AND. UR + snd .GE. 0.0D0) THEN
!C           3-centered rarefaction wave
            RHOS = RHOR * DEXP(-(UR + snd) / snd)
            US(I) = -snd
            VS(I) = VR
         ELSE
!C           Right state
            RHOS = RHOR
            US(I) = UR
            VS(I) = VR
         END IF

!C        Fluctuations
         AMDQ(I,1)  = F1(RHOS,US(I),VS(I)) - F1(RHOL,UL,VL)
         AMDQ(I,IU) = F2(RHOS,US(I),VS(I)) - F2(RHOL,UL,VL)
         AMDQ(I,IV) = F3(RHOS,US(I),VS(I)) - F3(RHOL,UL,VL)
         APDQ(I,1)  = F1(RHOR,UR,VR) - F1(RHOS,US(I),VS(I))
         APDQ(I,IU) = F2(RHOR,UR,VR) - F2(RHOS,US(I),VS(I))
         APDQ(I,IV) = F3(RHOR,UR,VR) - F3(RHOS,US(I),VS(I))
   enddo

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modified from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine NEWTON(X0,EPSLON,ROOT,ERR,RHOL,UL,VL,RHOR,UR,VR)
use common_params
implicit none
double precision::X0,EPSLON,ROOT
logical:: ERR
double precision:: EQ4QM,EQ4QMP
double precision:: RHOL,UL,VL,RHOR,UR,VR

integer,parameter::ITRMAX=10
integer:: ITRNUM
double precision:: XNEW,XOLD
 
      XOLD = X0
      DO ITRNUM = 1, ITRMAX
         XNEW = XOLD-EQ4QM(XOLD,RHOL,UL,VL,RHOR,UR,VR)/EQ4QMP(XOLD,RHOL,UL,VL,RHOR,UR,VR)

         IF (DABS((XNEW - XOLD) / XOLD) .LT. EPSLON) THEN
            ROOT = XNEW
            ERR = .FALSE.
            RETURN
         END IF
         XOLD = XNEW
      enddo   
      ROOT = XNEW
      write(*,*) root,RHOL,UL,VL,RHOR,UR,VR
      ERR = .TRUE.


end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Copy from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION EQ4QM(RHO,RHOL,UL,VL,RHOR,UR,VR)
IMPLICIT NONE
DOUBLE PRECISION RHO
DOUBLE PRECISION PHIL, PHIR
double precision RHOL,UL,VL
double precision RHOR,UR,VR

      EQ4QM = PHIL(RHO,RHOL,UL,VL) - PHIR(RHO,RHOR,UR,VR)

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modified from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION PHIL (RHO,RHOL,UL,VL)
      use common_params
      IMPLICIT NONE
      DOUBLE PRECISION RHO
      DOUBLE PRECISION TOL
      DOUBLE PRECISION EPS
      DOUBLE PRECISION RHOL, UL, VL

      TOL=1.d-6

      EPS = (RHO - RHOL) / RHOL
      IF (RHO .LE. RHOL) THEN
         IF (DABS(EPS) .GE. TOL) THEN
            PHIL = UL + snd * DLOG(RHOL / RHO)
         ELSE
            PHIL = UL - snd * EPS * (1.0D0 - 0.5D0 * EPS &
                                        + EPS * EPS / 3.0D0)
         END IF
      ELSE
         IF (DABS(EPS) .GE. TOL) THEN
            PHIL = UL - snd * (RHO - RHOL) / DSQRT(RHOL * RHO)
         ELSE
            PHIL = UL - snd * EPS * (1.0D0 - 0.5D0 * EPS &
                                        + 0.375D0 * EPS * EPS)
         END IF
      END IF
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modified from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION PHIR (RHO,RHOR,UR,VR)
      use common_params
      implicit none
      double precision:: RHO,RHOR,UR,VR
      DOUBLE PRECISION TOL
      DOUBLE PRECISION EPS

      TOL=1.d-6

      EPS = (RHO - RHOR) / RHOR
      IF (RHO .LE. RHOR) THEN
         IF (DABS(EPS) .GE. TOL) THEN
            PHIR = UR - snd * DLOG(RHOR / RHO)
         ELSE
            PHIR = UR + snd * EPS * (1.0D0 - 0.5D0 * EPS &
                                        + EPS * EPS / 3.0D0)
         END IF
      ELSE
         IF (DABS(EPS) .GE. TOL) THEN
            PHIR = UR + snd * (RHO - RHOR) / DSQRT(RHOR * RHO)
         ELSE
            PHIR = UR + snd * EPS * (1.0D0 - 0.5D0 * EPS &
                                        + 0.375D0 * EPS * EPS)
         END IF
      END IF
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modifed from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION EQ4QMP (RHO,RHOL,UL,VL,RHOR,UR,VR)
      use common_params
      IMPLICIT NONE
      DOUBLE PRECISION RHO
      double precision::RHOL,UL,VL
      double precision::RHOR,UR,VR
      DOUBLE PRECISION PHILP, PHIRP

      EQ4QMP = PHILP(RHO,RHOL,UL,VL) - PHIRP(RHO,RHOR,UR,VR)
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modified from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION PHILP (RHO,RHOL,UL,VL)
use common_params
      IMPLICIT NONE
      DOUBLE PRECISION RHO
      DOUBLE PRECISION TOL
      DOUBLE PRECISION EPS
      double precision::RHOL,UL,VL

      TOL=1.d-6

      EPS = (RHO - RHOL) / RHOL
      IF (RHO .LE. RHOL) THEN
         IF (DABS(EPS) .GE. TOL) THEN
            PHILP = -snd / RHO
         ELSE
            PHILP = -(snd / RHOL) * (1.0D0 - EPS + EPS * EPS)
         END IF
      ELSE
         IF (DABS(EPS) .GE. TOL) THEN
            PHILP = -0.5D0 * snd * (RHOL/RHO + 1.0D0) / DSQRT(RHOL * RHO)
         ELSE
            PHILP = -(snd / RHOL) * (1.0D0 - EPS + 1.125D0 * EPS * EPS)
         END IF
      END IF

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modified from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION PHIRP (RHO,RHOR,UR,VR)
use common_params
      IMPLICIT NONE
      DOUBLE PRECISION RHO
      double precision:: RHOR,UR,VR
      DOUBLE PRECISION TOL
      DOUBLE PRECISION EPS

      TOL=1.d-6

      EPS = (RHO - RHOR) / RHOR
      IF (RHO .LE. RHOR) THEN
         IF (DABS(EPS) .GE. TOL) THEN
            PHIRP = snd / RHO
         ELSE
            PHIRP = (snd / RHOR) * (1.0D0 - EPS + EPS * EPS)
         END IF
      ELSE
         IF (DABS(EPS) .GE. TOL) THEN
            PHIRP = 0.5D0 * snd * (RHOR / RHO + 1.0D0) / DSQRT(RHOR * RHO)
         ELSE
            PHIRP = (snd / RHOR) * (1.0D0 - EPS + 1.125D0 * EPS * EPS)
         END IF
      END IF
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modified from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LIMIT2(W,NX1D,NBC,NW,S)
use common_params
implicit none
integer::NX1D,NBC,NW
double precision::W(2-NBC:NX1D+NBC,NVAR,NW),S(2-NBC:NX1D+NBC,NW)
integer::I,K,L
double precision::DOTL,DOTR,NORM2,LIMIT
double precision::PHILIM

!C Subroutine Body
      DO L = 1, NW
         DOTR = 0.0D0
         DO  K = 1, NVAR
            DOTR = DOTR + W(0,K,L) * W(1,K,L)
         enddo

         DO I = 1, NX1D + 1
!C           Calculate the dot product and the squared norm of the waves.
            DOTL = DOTR
            DOTR = 0.0D0
            NORM2 = 0.0D0
            DO K = 1, NVAR
               DOTR = DOTR + W(I,K,L) * W(I+1,K,L)
               NORM2 = NORM2 + W(I,K,L) ** 2
            enddo

            IF (NORM2 .NE. 0.0D0) THEN
!C              Choose the upwind direction.
               IF (S(I,L) .GT. 0.0D0) THEN
                  LIMIT = PHILIM(DOTL / NORM2)
               ELSE
                  LIMIT = PHILIM(DOTR / NORM2)
               END IF
!C              Apply the limiter.
               DO K = 1, NVAR
                  W(I,K,L) = LIMIT * W(I,K,L)
               enddo
            END IF
        enddo 
      enddo   

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modified from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION PHILIM (THETA)
use common_params
implicit none
double precision::THETA
integer,parameter:: IL=5

!C Function Body
      IF (IL .EQ. 1) THEN
!C        Upwind
         PHILIM = 0.0D0
      ELSE IF (IL .EQ. 2) THEN
!C        Lax-Wendroff
         PHILIM = 1.0D0
      ELSE IF (IL .EQ. 3) THEN
!C        Beam-Warming
         PHILIM = THETA
      ELSE IF (IL .EQ. 4) THEN
!C        Fromm
         PHILIM = 0.5D0 * (1.0D0 + THETA)
      ELSE IF (IL .EQ. 5) THEN
!C        Minmod
         IF (1.0D0 .LE. THETA) THEN
            PHILIM = 1.0D0
         ELSE IF (0.0D0 .LT. THETA .AND. THETA .LT. 1.0D0) THEN
            PHILIM = THETA
         ELSE
            PHILIM = 0.0D0
         END IF
      ELSE IF (IL .EQ. 6) THEN
!C        Superbee
         PHILIM = DMAX1(0.0D0, DMIN1(1.0D0, 2.0D0 * THETA), &
                           DMIN1(2.0D0, THETA))
      ELSE IF (IL .EQ. 7) THEN
!C        MC
         PHILIM = DMAX1(0.0D0, DMIN1(0.5D0 * (1.0D0 + THETA), &
                                 2.0D0, 2.0D0 * THETA))
      ELSE IF (IL .EQ. 8) THEN
!C        van Leer
         PHILIM = (THETA + DABS(THETA)) / (1.0D0 + DABS(THETA))
      ELSE
!C        Error Section
         PRINT *, 'Unknown choice for the flux limiter.'
         STOP
      END IF

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! modified from Chao-Chin Yang's code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RPT2(NX1D,NBC,IXY,IMP,Q1D,ASDQ,BMASDQ,BPASDQ,US,VS)
use common_params
implicit none
integer::NX1D,NBC,IXY,IMP
double precision::Q1D(1-NBC:NX1D+NBC,NVAR)
double precision::ASDQ(2-NBC:NX1D+NBC,NVAR)
double precision::BMASDQ(2-NBC:NX1D+NBC,NVAR)
double precision::BPASDQ(2-NBC:NX1D+NBC,NVAR)
integer:: IU,IV
double precision::BETA(3),LAMBDA(3)
integer::I,K
double precision::CP(3),CM(3),TEMP
double precision::US(2-NBC:NX1D+NBC),VS(2-NBC:NX1D+NBC)


!C Check the normal direction.
      IF (IXY .EQ. 1) THEN
         IU = 2
         IV = 3
      ELSE IF (IXY .EQ. 2) THEN
         IU = 3
         IV = 2
      ELSE
         PRINT *, 'Unknown sweep code.'
         STOP
      END IF

!C Decompose ASDQ at each interface.
      TEMP = 2.0D0 * snd
      DO I = 2 - NBC, NX1D + NBC
         LAMBDA(1) = VS(I) - snd
         LAMBDA(2) = VS(I)
         LAMBDA(3) = VS(I) + snd
         BETA(1) = (LAMBDA(3) * ASDQ(I,1) - ASDQ(I,IV)) / TEMP
         BETA(2) = -US(I) * ASDQ(I,1) + ASDQ(I,IU)
         BETA(3) = (-LAMBDA(1) * ASDQ(I,1) + ASDQ(I,IV)) / TEMP
         DO K = 1, 3
            CP(K) = DMAX1(LAMBDA(K), 0.0D0) * BETA(K)
            CM(K) = DMIN1(LAMBDA(K), 0.0D0) * BETA(K)
         enddo
         BPASDQ(I,1) = CP(1) + CP(3)
         BPASDQ(I,IU) = US(I) * BPASDQ(I,1) + CP(2)
         BPASDQ(I,IV) = LAMBDA(1) * CP(1) + LAMBDA(3) * CP(3)
         BMASDQ(I,1) = CM(1) + CM(3)
         BMASDQ(I,IU) = US(I) * BMASDQ(I,1) + CM(2)
         BMASDQ(I,IV) = LAMBDA(1) * CM(1) + LAMBDA(3) * CM(3)
      enddo

end subroutine
