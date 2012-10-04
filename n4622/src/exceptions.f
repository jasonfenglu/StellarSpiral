      SUBROUTINE newton_exception(i,j,d)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER                           ::i,j
      character                         ::d
      INCLUDE 'mpif.h'
      INCLUDE '../parameters'
      call vars(xmin,xmax,ymin,ymax,
     &          dcfl,ncir,psrl,ichk,
     &          fstr,fstp,nstp,unk1)
      write(*,*)'myid is ',myid
      if (d .eq. 'x')then
              write(*,*)'grid (i,j) = ',i,j
              write(*,*)'position: (x,y) = ',x(i),y(j)
      elseif (d .eq. 'y')then
              write(*,*)'grid (i,j) = ',j,i
              write(*,*)'position: (x,y) = ',x(j),y(i)
      endif
       
      STOP
     
      ENDSUBROUTINE newton_exception
