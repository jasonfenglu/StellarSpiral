!$Id: pgplot_2dmesh.f 156 2011-08-05 10:46:38Z ccfeng $
        MODULE MATPLOT
        CONTAINS

      SUBROUTINE drawmat(F,MXI,MXJ)
      !DEC$ ATTRIBUTES ALIAS : 'drawmat' :: drawmat
C-----------------------------------------------------------------------
C Test program for PGPLOT: test of imaging routine PGIMAG and associated
C routines PGWEDG and PGCTAB.
C-----------------------------------------------------------------------
      INTEGER PGOPEN
      INTEGER   MXI, MXJ
c     PARAMETER (MXI=64, MXJ=64)
      INTEGER I, L, C1, C2, NC
      REAL                              :: F(MXI,MXJ)
      REAL FMIN,FMAX,TR(6), CONTRA, BRIGHT, ANGLE, C, S, ALEV(1)
      CHARACTER*16 VAL

c     CALL print_matrix(F,MXJ,MXI)


C
C Introduction.
C
c     WRITE(*,*)'Demonstration of PGIMAG and associated routines.'
c     WRITE(*,*)'This program requires a device with color capability.'
c     WRITE(*,*)'On an interactive device, you can modify the color map'
c     WRITE(*,*)'used for the image.'
c     WRITE(*,*)
C
C Open device for graphics.
C
!     IF (PGOPEN('?') .LT. 1) STOP
!     CALL PGQINF('TYPE', VAL, L)
!     WRITE (*,*) 'PGPLOT device type: ', VAL(1:L)
!     CALL PGQCIR(C1, C2)
!     NC = MAX(0, C2-C1+1)
!     WRITE (*,*) 'Number of color indices used for image: ', NC
!     IF (NC .LT.8) THEN 
!        WRITE (*,*) 'Not enough colors available on this device'
!        STOP
!     ELSE
!        WRITE (*,*)
!     END IF
C
C Compute a suitable function in array F.
C

        

      i = PGOPEN('/XWINDOW')
c     i = PGOPEN('/GIF')
c     i = PGOPEN('?')
c     CALL FUNC(F, MXI, MXJ, FMIN, FMAX)
C
C-----------------------------------------------------------------------
C Example 1: simple transformation matrix
C-----------------------------------------------------------------------
C
C Set the coordinate transformation matrix: 
C world coordinate = pixel number.
C
      TR(1) = 0.0
      TR(2) = 1.0
      TR(3) = 0.0
c     TR(4) = 0.0
      TR(4) = 0.0
      TR(5) = 0.0
      TR(6) = 1.0
C
C Clear the screen. Set up window and viewport.
C
      CALL PGPAGE
      CALL SETVP
      CALL PGWNAD(0.0, 1.0+MXI, 0.0, 1.0+MXJ)
C
C Set up the color map.
C
      BRIGHT = 0.5
      CONTRA  = 1.0
      CALL PALETT(1, CONTRA, BRIGHT)
C
C Draw the map with PGIMAG.  
C
      FMIN = MINVAL(F)
      FMAX = MAXVAL(F)
      WRITE(*,*)'plot min max',FMIN,FMAX
      !FMAX = 4.
      !FMIN = -4.
      CALL PGIMAG(F,MXI,MXJ,1,MXI,1,MXJ,FMIN,FMAX,TR)
c     CALL PGIMAG(F,MXI,MXJ,1,MXI,1,MXJ,0.,255.,TR)

C
C Annotate the plot.
C
      CALL PGMTXT('t',1.0,0.0,0.0,'PGIMAG, PGWEDG, and PGCTAB')
      CALL PGSCH(0.6)
      CALL PGBOX('bcntsi',0.0,0,'bcntsiv',0.0,0)
      CALL PGMTXT('b',3.0,1.0,1.0,'pixel number')
C
C Draw a wedge.
C
      CALL PGWEDG('BI', 4.0, 5.0, FMIN, FMAX, 'pixel value')
      CALL PGSCH(1.0)
C
C If the device has a cursor, allow user to fiddle with color table.
C
c     CALL PGQINF('CURSOR', VAL, L)
c     IF (VAL(:L).EQ.'YES') THEN
c        CALL FIDDLE
c        CALL PGASK(.FALSE.)
c     END IF
C
C Close the device and exit.
C
      CALL PGEND
C-----------------------------------------------------------------------
      ENDSUBROUTINE drawmat

      SUBROUTINE drawvec(Fx,Fy,MXI,MXJ)
      !DEC$ ATTRIBUTES ALIAS : 'drawmat' :: drawmat
C-----------------------------------------------------------------------
C Test program for PGPLOT: test of imaging routine PGIMAG and associated
C routines PGWEDG and PGCTAB.
C-----------------------------------------------------------------------
      INTEGER PGOPEN
      INTEGER   MXI, MXJ
c     PARAMETER (MXI=64, MXJ=64)
      INTEGER I, L, C1, C2, NC
      REAL                              :: Fx(MXI,MXJ),Fy(MXI,MXJ)
      REAL FMIN,FMAX,TR(6), CONTRA, BRIGHT, ANGLE, C, S, ALEV(1)
      CHARACTER*16 VAL

c     CALL print_matrix(F,MXJ,MXI)


C
C Introduction.
C
c     WRITE(*,*)'Demonstration of PGIMAG and associated routines.'
c     WRITE(*,*)'This program requires a device with color capability.'
c     WRITE(*,*)'On an interactive device, you can modify the color map'
c     WRITE(*,*)'used for the image.'
c     WRITE(*,*)
C
C Open device for graphics.
C
!     IF (PGOPEN('?') .LT. 1) STOP
!     CALL PGQINF('TYPE', VAL, L)
!     WRITE (*,*) 'PGPLOT device type: ', VAL(1:L)
!     CALL PGQCIR(C1, C2)
!     NC = MAX(0, C2-C1+1)
!     WRITE (*,*) 'Number of color indices used for image: ', NC
!     IF (NC .LT.8) THEN 
!        WRITE (*,*) 'Not enough colors available on this device'
!        STOP
!     ELSE
!        WRITE (*,*)
!     END IF
C
C Compute a suitable function in array F.
C

        

c     i = PGOPEN('/XWINDOW')
c     i = PGOPEN('/GIF')
      i = PGOPEN('?')
c     CALL FUNC(F, MXI, MXJ, FMIN, FMAX)
C
C-----------------------------------------------------------------------
C Example 1: simple transformation matrix
C-----------------------------------------------------------------------
C
C Set the coordinate transformation matrix: 
C world coordinate = pixel number.
C
      TR(1) = 0.0
      TR(2) = 1.0
      TR(3) = 0.0
      TR(4) = 0.0
      TR(5) = 0.0
      TR(6) = 1.0
C
C Clear the screen. Set up window and viewport.
C
      CALL PGPAGE
      CALL SETVP
      CALL PGWNAD(0.0, 1.0+MXI, 0.0, 1.0+MXJ)
C
C Set up the color map.
C
      BRIGHT = 0.5
      CONTRA  = 5.0
      CALL PALETT(1, CONTRA, BRIGHT)
C
C Draw the map with PGIMAG.  
C
      FMIN = MINVAL(Fx)
      FMAX = MAXVAL(Fx)
      WRITE(*,*)'plot min max',FMIN,FMAX
      !FMAX = 4.
      !FMIN = -4.
      !CALL PGIMAG(F,MXI,MXJ,1,MXI,1,MXJ,FMIN,FMAX,TR)
      CALL PGSCH(0.6)
      CALL PGVECT(Fx,Fy,MXI,MXJ,1,MXI,1,MXJ,0.0,0,TR,0.)
c     CALL PGIMAG(F,MXI,MXJ,1,MXI,1,MXJ,0.,255.,TR)

C
C Annotate the plot.
C
      CALL PGMTXT('t',1.0,0.0,0.0,'PGIMAG, PGWEDG, and PGCTAB')
      CALL PGSCH(0.6)
      CALL PGBOX('bcntsi',0.0,0,'bcntsiv',0.0,0)
      CALL PGMTXT('b',3.0,1.0,1.0,'pixel number')
C
C Draw a wedge.
C
      CALL PGWEDG('BI', 4.0, 5.0, FMIN, FMAX, 'pixel value')
      CALL PGSCH(1.0)
C
C If the device has a cursor, allow user to fiddle with color table.
C
c     CALL PGQINF('CURSOR', VAL, L)
c     IF (VAL(:L).EQ.'YES') THEN
c        CALL FIDDLE
c        CALL PGASK(.FALSE.)
c     END IF
C
C Close the device and exit.
C
      CALL PGEND
C-----------------------------------------------------------------------
      ENDSUBROUTINE drawvec

      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
C-----------------------------------------------------------------------
C Set a "palette" of colors in the range of color indices used by
C PGIMAG.
C-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      ENDSUBROUTINE PALETT

      SUBROUTINE SETVP
C-----------------------------------------------------------------------
C Set the viewport, allowing margins around the edge for annotation.
C (This is similar in effect to PGVSTD, but has different margins.)
C The routine determines the view-surface size and allocates margins
C as fractions of the minimum of width and height.
C-----------------------------------------------------------------------
      REAL D, VPX1, VPX2, VPY1, VPY2
C
      CALL PGSVP(0.0, 1.0, 0.0, 1.0)
      CALL PGQVP(1, VPX1, VPX2, VPY1, VPY2)
      D = MIN(VPX2-VPX1, VPY2-VPY1)/40.0
      VPX1 = VPX1 + 5.0*D
      VPX2 = VPX2 - 2.0*D
      VPY1 = VPY1 + 8.0*D
      VPY2 = VPY2 - 2.0*D
      CALL PGVSIZ(VPX1, VPX2, VPY1, VPY2)
      ENDSUBROUTINE SETVP

      SUBROUTINE FIDDLE
C
      INTEGER P, IER, PGCURS
      REAL CONTRA, BRIGHT, X, Y, SIGN
      REAL X1, Y1, X2, Y2, B1, B2, C1, C2
      CHARACTER CH
C
      WRITE (*,*) 'Use cursor to adjust color table:'
      WRITE (*,*) ' Keys 1,2,3,4,5 select different palettes'
      WRITE (*,*) ' Key P cycles through available palettes'
      WRITE (*,*) ' Key F adjusts contrast and brightness, with'
      WRITE (*,*) '  cursor x position setting brightness [0.0 - 1.0]'
      WRITE (*,*) '   and y position setting contrast [0.0 - 10.0]'
      WRITE (*,*) '  (Hold down F key while moving cursor to change'
      WRITE (*,*) '  contrast and brightness continuously)'
      WRITE (*,*) ' Key C resets contrast=1.0, brightness=0.5'
      WRITE (*,*) ' Key - reverses color palette'
      WRITE (*,*) ' Key X or right mouse button exits program' 
C
      P = 2
      CONTRA = 1.0
      BRIGHT = 0.5
      X = 0.5
      Y = 1.0
      SIGN = +1.0
C
      CALL PGQWIN(X1, X2, Y1, Y2)
      B1 = 0.0
      B2 = 1.0
      C1 = 0.0
      C2 = 10.0
      CALL PGSWIN(B1, B2, C1, C2)
 10   IER = PGCURS(X, Y, CH)
      IF (CH.EQ.CHAR(0) .OR. CH.EQ.'x' .OR. CH.EQ.'X') THEN
         CALL PGSWIN(X1, X2, Y1, Y2)
         RETURN
      ELSE IF (CH.EQ.'F' .OR. CH.EQ.'f') THEN
         BRIGHT = MAX(B1, MIN(B2,X))
         CONTRA = MAX(C1, MIN(C2,Y))
      ELSE IF (CH.EQ.'C' .OR. CH.EQ.'c') THEN
         CONTRA = 1.0
         Y = 1.0
         BRIGHT = 0.5
         X = 0.5
      ELSE IF (CH.EQ.'-') THEN
         SIGN = -SIGN
      ELSE IF (CH.EQ.'1') THEN
         P = 1
      ELSE IF (CH.EQ.'2') THEN
         P = 2
      ELSE IF (CH.EQ.'3') THEN
         P = 3
      ELSE IF (CH.EQ.'4') THEN
         P = 4
      ELSE IF (CH.EQ.'5') THEN
         P = 5
      ELSE IF (CH.EQ.'P' .OR. CH.EQ.'p') THEN
         P = 1 + MOD(P,5)
      END IF
      CALL PALETT(P, SIGN*CONTRA, BRIGHT)
      GOTO 10
      ENDSUBROUTINE FIDDLE

      ENDMODULE MATPLOT
