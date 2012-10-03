      PROGRAM PGDEM3
C-----------------------------------------------------------------------
C Demonstration program for PGPLOT contouring routines.
C-----------------------------------------------------------------------
      INTEGER PGBEG
      WRITE (*,'(A)') ' Demonstration of PGPLOT contouring routines'
C
C Call PGBEG to initiate PGPLOT and open the output device; PGBEG
C will prompt the user to supply the device name and type.
C
      IF (PGBEG(0,'?',1,1) .NE. 1) STOP
C
C Call the demonstration subroutines.
C
      WRITE (*,'(A)') ' Routine PGCONF'
      CALL PGEXX1
C
C Finally, call PGEND to terminate things properly.
C
      CALL PGEND
C-----------------------------------------------------------------------
      END


      SUBROUTINE PGEXX1
C-----------------------------------------------------------------------
C Demonstration of contouring routine PGCONF.
C-----------------------------------------------------------------------
      INTEGER NX, NY, NC
      PARAMETER (NX=51, NY=51, NC=9)
      INTEGER I, J
      REAL Z(NX,NY),TR(6), R
      REAL X, Y, XMIN, XMAX, YMIN, YMAX, DX, DY, MU, C(NC)
      DATA C /3.0, 3.2, 3.5, 3.6, 3.766413, 4.0 ,5.0, 10.0, 100.0/     
C
C Compute a suitable function. This is the test function used by
C W. V. Snyder, Algorithm 531, Contour Plotting, ACM Trans. Math.
C Softw. v.4, pp.290-294 (1978).
C
      XMIN = -2.0
      XMAX = 2.0
      YMIN =-2.0
      YMAX = 2.0
      MU = 0.3
      DX = (XMAX-XMIN)/FLOAT(NX-1)                                      
      DY = (YMAX-YMIN)/FLOAT(NY-1)
      TR(1) = XMIN - DX
      TR(2) = DX
      TR(3) = 0.0
      TR(4) = YMIN - DY
      TR(5) = 0.0
      TR(6) = DY
      DO 20 I=1,NX
         X = TR(1) + I*TR(2)
         DO 10 J=1,NY     
            Y = TR(4) + J*TR(6)
            Z(I,J) = (1.0-MU)*(2.0/SQRT((X-MU)**2+Y**2)+(X-MU)**2+Y**2)   
     *           + MU*(2.0/SQRT((X+1.0-MU)**2+Y**2)+(X+1.0-MU)**2+Y**2)      
 10      CONTINUE                                   
   20 CONTINUE                                                          
C
C Clear the screen. Set up window and viewport.
C
      CALL PGPAGE
      CALL PGVSTD(0.05,0.95,0.05,0.95)
      CALL PGWNAD(XMIN, XMAX, YMIN, YMAX)
C
C Fill contours with PGCONF.
C
c     CALL PGSFS(1)
c     DO 30 I=1, NC-1
c        R = 0.5+0.5*REAL(I-1)/REAL(NC-1)
c        CALL PGSCR(I+10, R, R, R)
c        CALL PGSCI(I+10)
c        CALL PGCONF(Z,NX,NY,1,NX,1,NY,C(I),C(I+1),TR)
c30   CONTINUE
      CALL PGIMAG(Z,NX,NY,1,NX,1,NY,C(1),C(NC),TR)
C
C Draw the contour lines with PGCONT.
C
      CALL PGSCI(3)
c     CALL PGCONT(Z,NX,NY,1,NX,1,NY,C,NC,TR)
C
C Labels and box.
C
      CALL PGSCI(1)
      CALL PGSCH(0.6)
      CALL PGBOX('bctsin',1.0,10,'bctsinv',1.0,10)
      CALL PGSCH(1.0)
      CALL PGMTXT('t',1.0,0.0,0.0,'Contour filling using PGCONF')
C
      END
