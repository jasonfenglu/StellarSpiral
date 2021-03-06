*DECK DPPGQ8
      SUBROUTINE DPPGQ8 (FUN, LDC, C, XI, LXI, KK, ID, A, B, INPPV, ERR,
     +   ANS, IERR)
C***BEGIN PROLOGUE  DPPGQ8
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DPFQAD
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (PPGQ8-S, DPPGQ8-D)
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract    **** A DOUBLE PRECISION routine ****
C
C        DPPGQ8, a modification of GAUS8, integrates the
C        product of FUN(X) by the ID-th derivative of a spline
C        DPPVAL(LDC,C,XI,LXI,KK,ID,X,INPPV)  between limits A and B.
C
C     Description of Arguments
C
C      Input-- FUN,C,XI,A,B,ERR are DOUBLE PRECISION
C        FUN - Name of external function of one argument which
C              multiplies DPPVAL.
C        LDC - Leading dimension of matrix C, LDC .GE. KK
C        C   - Matrix of Taylor derivatives of dimension at least
C              (K,LXI)
C        XI  - Breakpoint vector of length LXI+1
C        LXI - Number of polynomial pieces
C        KK  - Order of the spline, KK .GE. 1
C        ID  - Order of the spline derivative, 0 .LE. ID .LE. KK-1
C        A   - Lower limit of integral
C        B   - Upper limit of integral (may be less than A)
C        INPPV- Initialization parameter for DPPVAL
C        ERR - Is a requested pseudorelative error tolerance.  Normally
C              pick a value of ABS(ERR) .LT. 1D-3.  ANS will normally
C              have no more error than ABS(ERR) times the integral of
C              the absolute value of FUN(X)*DPPVAL(LDC,C,XI,LXI,KK,ID,X,
C              INPPV).
C
C
C      Output-- ERR,ANS are DOUBLE PRECISION
C        ERR - Will be an estimate of the absolute error in ANS if the
C              input value of ERR was negative.  (ERR Is unchanged if
C              the input value of ERR was nonnegative.)  The estimated
C              error is solely for information to the user and should
C              not be used as a correction to the computed integral.
C        ANS - Computed value of integral
C        IERR- A status code
C            --Normal Codes
C               1 ANS most likely meets requested error tolerance,
C                 or A=B.
C              -1 A and B are too nearly equal to allow normal
C                 integration.  ANS is set to zero.
C            --Abnormal Code
C               2 ANS probably does not meet requested error tolerance.
C
C***SEE ALSO  DPFQAD
C***ROUTINES CALLED  D1MACH, DPPVAL, I1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
C***END PROLOGUE  DPPGQ8
C
      INTEGER ID,IERR,INPPV,K,KK,KML,KMX,L,LDC,LMN,LMX,LR,LXI,MXL,
     1 NBITS, NIB, NLMN, NLMX
      INTEGER I1MACH
      DOUBLE PRECISION A,AA,AE,ANIB,ANS,AREA,B,BE,C,CC,EE,EF,EPS,ERR,
     1 EST,GL,GLR,GR,HH,SQ2,TOL,VL,VR,W1, W2, W3, W4, XI, X1,
     2 X2, X3, X4, X, H
      DOUBLE PRECISION D1MACH, DPPVAL, G8, FUN
      DIMENSION XI(*), C(LDC,*)
      DIMENSION AA(60), HH(60), LR(60), VL(60), GR(60)
      SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2, NLMN, KMX, KML
      DATA X1, X2, X3, X4/
     1     1.83434642495649805D-01,     5.25532409916328986D-01,
     2     7.96666477413626740D-01,     9.60289856497536232D-01/
      DATA W1, W2, W3, W4/
     1     3.62683783378361983D-01,     3.13706645877887287D-01,
     2     2.22381034453374471D-01,     1.01228536290376259D-01/
      DATA SQ2/1.41421356D0/
      DATA NLMN/1/,KMX/5000/,KML/6/
      G8(X,H)=
     1   H*((W1*(FUN(X-X1*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X-X1*H,INPPV)
     2          +FUN(X+X1*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X+X1*H,INPPV))
     3      +W2*(FUN(X-X2*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X-X2*H,INPPV)
     4          +FUN(X+X2*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X+X2*H,INPPV)))
     5     +(W3*(FUN(X-X3*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X-X3*H,INPPV)
     6          +FUN(X+X3*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X+X3*H,INPPV))
     7      +W4*(FUN(X-X4*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X-X4*H,INPPV)
     8          +FUN(X+X4*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X+X4*H,INPPV))))
C
C     INITIALIZE
C
C***FIRST EXECUTABLE STATEMENT  DPPGQ8
      K = I1MACH(14)
      ANIB = D1MACH(5)*K/0.30102000D0
      NBITS = INT(ANIB)
      NLMX = MIN((NBITS*5)/8,60)
      ANS = 0.0D0
      IERR = 1
      BE = 0.0D0
      IF (A.EQ.B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B.EQ.0.0D0) GO TO 10
      IF (SIGN(1.0D0,B)*A.LE.0.0D0) GO TO 10
      CC = ABS(1.0D0-A/B)
      IF (CC.GT.0.1D0) GO TO 10
      IF (CC.LE.0.0D0) GO TO 140
      ANIB = 0.5D0 - LOG(CC)/0.69314718D0
      NIB = INT(ANIB)
      LMX = MIN(NLMX,NBITS-NIB-7)
      IF (LMX.LT.1) GO TO 130
      LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0D0**(5-NBITS))/2.0D0
      IF (ERR.EQ.0.0D0) TOL = SQRT(D1MACH(4))
      EPS = TOL
      HH(1) = (B-A)/4.0D0
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0D0*HH(L),2.0D0*HH(L))
      K = 8
      AREA = ABS(EST)
      EF = 0.5D0
      MXL = 0
C
C     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
C
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0D0*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
      GLR = GL + GR(L)
      EE = ABS(EST-GLR)*EF
      AE = MAX(EPS*AREA,TOL*ABS(GLR))
      IF (EE-AE) 40, 40, 50
   30 MXL = 1
   40 BE = BE + (EST-GLR)
      IF (LR(L)) 60, 60, 80
C
C     CONSIDER THE LEFT HALF OF THIS LEVEL
C
   50 IF (K.GT.KMX) LMX = KML
      IF (L.GE.LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5D0
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5D0
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
C
C     PROCEED TO RIGHT HALF AT THIS LEVEL
C
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0D0*HH(L)
      GO TO 20
C
C     RETURN ONE LEVEL
C
   80 VR = GLR
   90 IF (L.LE.1) GO TO 120
      L = L - 1
      EPS = EPS*2.0D0
      EF = EF*SQ2
      IF (LR(L)) 100, 100, 110
  100 VL(L) = VL(L+1) + VR
      GO TO 70
  110 VR = VL(L+1) + VR
      GO TO 90
C
C      EXIT
C
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (ABS(BE).LE.2.0D0*TOL*AREA)) GO TO 140
      IERR = 2
      CALL XERMSG ('SLATEC', 'DPPGQ8',
     +   'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.', 3, 1)
      GO TO 140
  130 IERR = -1
      CALL XERMSG ('SLATEC', 'DPPGQ8',
     +   'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL ' //
     +   'INTEGRATION.  ANSWER IS SET TO ZERO, AND IERR=-1.', 1, -1)
  140 CONTINUE
      IF (ERR.LT.0.0D0) ERR = BE
      RETURN
      END
