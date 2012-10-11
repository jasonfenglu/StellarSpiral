      PROGRAM PGDEM3
      USE ImageTool
C-----------------------------------------------------------------------
C Demonstration program for PGPLOT contouring routines.
C-----------------------------------------------------------------------
      INTEGER PGBEG
      double precision, ALLOCATABLE            ::dat(:,:)
      INTEGER                      ::naxes(2)

      type (image) M81Image
C
C     FIND IMAGE SIZE AND COPY INTO ARRAY
C

      CALL readimagesize(naxes)
      ALLOCATE(dat(naxes(1),naxes(2)))
      CALL readimage(dat,naxes)
      
      CALL InitImage(M81Image,dat,naxes,(/1.d0,1.d0/),(/10.d0,10.d0/))
      DEALLOCATE(dat)

      write(*,*)'?????????'
!     write(*,*)'to real',PixelValue((/1.d0,1.d0/),M81Image)
!     write(*,*)'to real',PixelValue()

      
C
C Call PGBEG to initiate PGPLOT and open the output device; PGBEG
C will prompt the user to supply the device name and type.
C
      IF (PGBEG(0,'?',1,1) .NE. 1) STOP
C
C Call the demonstration subroutines.
C
      CALL plotcont(M81Image%imdata,M81Image%DataSize)
      CALL PGIDEN
C
C Finally, call PGEND to terminate things properly.
C
      CALL PGEND

      CALL DeallocateImage(M81Image)
        


      END

      subroutine readimage(dat,naxes)

C  Read a FITS image and determine the minimum and maximum pixel value.
C  Rather than reading the entire image in
C  at once (which could require a very large array), the image is read
C  in pieces, 100 pixels at a time.  

      integer status,unit,readwrite,blocksize,naxes(2),nfound
      integer group,firstpix,nbuffer,npixels,i
      real    dat(naxes(1)*naxes(2))
      real    datamin,datamax,nullval,buffer(100)
      logical anynull
      character filename*80


C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

C  Open the FITS file previously created by WRITEIMAGE
c     filename='ATESTFILEZ.FITS'
      filename='M81.fits'
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

C  Determine the size of the image.
      call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

C  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the NAXISn keywords.'
          return
       end if

C  Initialize variables
      npixels=naxes(1)*naxes(2)
      group=1
      firstpix=1
      nullval=-999
      datamin=1.0E30
      datamax=-1.0E30

      do while (npixels .gt. 0)
C         read up to 100 pixels at a time 
          nbuffer=min(100,npixels)
      
          call ftgpve(unit,group,firstpix,nbuffer,nullval,
     &            buffer,anynull,status)

          dat(firstpix:firstpix+nbuffer) = buffer

C         increment pointers and loop back to read the next group of pixels
          npixels=npixels-nbuffer
          firstpix=firstpix+nbuffer
      end do


C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

C  Check for any error, and if so print out error messages.
C  The PRINTERROR subroutine is listed near the end of this file.
      if (status .gt. 0)call printerror(status)
      end

      subroutine readimagesize(naxes)

C  Read a FITS image and determine the minimum and maximum pixel value.
C  Rather than reading the entire image in
C  at once (which could require a very large array), the image is read
C  in pieces, 100 pixels at a time.  

      integer status,unit,readwrite,blocksize,naxes(2),nfound
      integer group,firstpix,nbuffer,npixels,i
      real datamin,datamax,nullval,buffer(100)
      logical anynull
      character filename*80


C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

C  Open the FITS file previously created by WRITEIMAGE
      filename='M81.fits'
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

C  Determine the size of the image.
      call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

C  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the NAXISn keywords.'
          return
       end if

C  Initialize variables
      npixels=naxes(1)*naxes(2)

C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

C  Check for any error, and if so print out error messages.
C  The PRINTERROR subroutine is listed near the end of this file.
      if (status .gt. 0)call printerror(status)
      end

      SUBROUTINE plotcont(dat,naxes)
      INTEGER, INTENT(IN)       ::      naxes(2)
      INTEGER, INTENT(IN)       ::      dat(naxes(1),naxes(2))
      INTEGER                   ::      NX,NY
      INTEGER                   ::      I, J
      REAL                      ::      TR(6), R
      REAL X, Y, XMIN, XMAX, YMIN, YMAX, DX, DY, MU


      NX = naxes(1)
      NY = naxes(2)

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
      CALL PGVSTD(0.05,0.95,0.05,0.95)
      CALL PGWNAD(0.0, REAL(NX), 0.0, REAL(NY))

C
C Fill contours with PGGRAY.
C

        write(*,*)'max',maxval(dat),minval(dat)
        CALL PGGRAY(dat,NX,NY,1,NX,1,NY,4.0,-0.0,TR) 
C
C Draw the contour lines with PGCONT.
C
c     CALL PGSCI(3)
c     CALL PGCONT(Z,NX,NY,1,NX,1,NY,C,NC,TR)
C
C Labels and box.
C
c     CALL PGSCI(1)
c     CALL PGSCH(0.6)
c     CALL PGBOX('bctsin',0.0,10,'bctsinv',1.0,10)
c     CALL PGSCH(1.0)
c     CALL PGMTXT('t',1.0,0.0,0.0,'Contour filling using PGCONF')

      END

      subroutine printerror(status)

C  This subroutine prints out the descriptive text corresponding to the
C  error status value and prints out the contents of the internal
C  error message stack generated by FITSIO whenever an error occurs.

      integer status
      character errtext*30,errmessage*80

C  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return

C  The FTGERR subroutine returns a descriptive 30-character text string that
C  corresponds to the integer error status number.  A complete list of all
C  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

C  FITSIO usually generates an internal stack of error messages whenever
C  an error occurs.  These messages provide much more information on the
C  cause of the problem than can be provided by the single integer error
C  status value.  The FTGMSG subroutine retrieves the oldest message from
C  the stack and shifts any remaining messages on the stack down one
C  position.  FTGMSG is called repeatedly until a blank message is
C  returned, which indicates that the stack is empty.  Each error message
C  may be up to 80 characters in length.  Another subroutine, called
C  FTCMSG, is available to simply clear the whole error message stack in
C  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      end
