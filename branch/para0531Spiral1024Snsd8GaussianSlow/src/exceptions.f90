SUBROUTINE newton_exception(i,j,d)
use common_params
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
INTEGER                           ::i,j
character                         ::d
INCLUDE 'mpif.h'
write(0,*)'myid is ',myid
if (d .eq. 'x')then
        write(0,*)'grid (i,j) = ',i,j
        write(0,*)'position: (x,y) = ',x_loc(i),y_loc(j)
elseif (d .eq. 'y')then
        write(0,*)'grid (i,j) = ',j,i
        write(0,*)'position: (x,y) = ',x_loc(j),y_loc(i)
endif
 
STOP

ENDSUBROUTINE newton_exception
