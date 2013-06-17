MODULE MATH
IMPLICIT NONE
TYPE    typintplt2
        DOUBLE PRECISION,ALLOCATABLE    ::dat(:,:)
        DOUBLE PRECISION,ALLOCATABLE    ::x(:)
        DOUBLE PRECISION,ALLOCATABLE    ::y(:)
        INTEGER                         ::nx,ny
        DOUBLE PRECISION                ::dx,dy,minx,miny,maxx,maxy
        CONTAINS
        PROCEDURE                       ::init=>typintplt2_init
        PROCEDURE                       ::find=>typintplt2_find
        PROCEDURE                       ::dismiss=>typintplt2_deallocate
ENDTYPE
!!interpolation 1-d
INTERFACE intplt1
        MODULE PROCEDURE intplt
        MODULE PROCEDURE cintplt
ENDINTERFACE
CONTAINS

FUNCTION typintplt2_find(this,x,y)
IMPLICIT NONE
class(typintplt2)                       ::this
DOUBLE PRECISION,INTENT(IN)             ::x,y
DOUBLE PRECISION                        ::xx,yy
DOUBLE PRECISION                        ::typintplt2_find
DOUBLE PRECISION                        ::intb(4)
INTEGER                                 ::k=0,l=0

intb = 0.d0
IF(x<this.minx.or.x>this.maxx.or.y<this.miny.or.y>this.maxy)THEN
        typintplt2_find = 0.d0
        write(6,*)'[math]: Intepolation of 2D out of data range.'
        write(6,*)x,y
ELSE
        !! find left and bottom grid
        DO k = 1, size(this.x)
                if(this.x(k)>x)EXIT
        ENDDO
        DO l = 1, size(this.y)
                if(this.y(l)>y)EXIT
        ENDDO
        k = k - 1
        l = l - 1

        xx = (x - this.x(k))/this.dx
        yy = (y - this.y(l))/this.dy

        intb(1) = this.dat(k,l)
        intb(2) = this.dat(k+1,l) - intb(1)
        intb(3) = this.dat(k,l+1) - intb(1)
        intb(4) = -sum(intb(:)) + this.dat(k+1,l+1)
        typintplt2_find = intb(1) + intb(2)*xx + intb(3)*yy &
                        + intb(4)*xx*yy
ENDIF

ENDFUNCTION

SUBROUTINE typintplt2_deallocate(this)
IMPLICIT NONE
class(typintplt2)                       ::this
DEALLOCATE(this.dat)
DEALLOCATE(this.x)
DEALLOCATE(this.y)
ENDSUBROUTINE

SUBROUTINE typintplt2_init(this,dat,x,y)
IMPLICIT NONE
class(typintplt2)                       ::this
DOUBLE PRECISION                        ::dat(:,:)
DOUBLE PRECISION                        ::x(:)
DOUBLE PRECISION                        ::y(:)

this.nx = size(x)
this.ny = size(y)
IF(.not.ALLOCATED(this.dat))ALLOCATE(this.dat(this.nx,this.ny))
IF(.not.ALLOCATED(this.x))ALLOCATE(this.x(this.nx))
IF(.not.ALLOCATED(this.y))ALLOCATE(this.y(this.ny))
this.x = x
this.y = y
this.dat = dat
this.dx = x(2)-x(1)
this.dy = y(2)-y(1)
this.minx = minval(this.x)
this.miny = minval(this.y)
this.maxx = maxval(this.x)
this.maxy = maxval(this.y)
ENDSUBROUTINE

FUNCTION intplt(dat,rs,r)
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)     ::dat(:),rs(:)
DOUBLE PRECISION                ::intplt
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE PRECISION                ::X(4),Y(4),C(4)
INTEGER,PARAMETER               ::N = 4
INTEGER                         ::NDER = 0,IERR
DOUBLE PRECISION                ::XX,YFIT,YP,WORK
INTEGER                         ::i

!If r is inside first four sample points
if(rs(4).gt.r)then
        Y(1:4) = dat(1:4)
        X(1:4) = rs(1:4)
!       write(0,*)'less than 4!!!!!!'
!       write(0,*)Y
!       write(0,*)'!!!'
elseif(r.lt.maxval(rs))then
        DO i = 4, size(dat)
                if(rs(i).gt.r)then
                        exit
                endif
        ENDDO
        Y(1:4) = dat(i-2:i+1)
        X(1:4) = rs(i-2:i+1)
ELSE
        i = maxloc(rs,1)
        Y(1:4) = dat(i-3:i)
        X(1:4) =  rs(i-3:i)
ENDIF        

CALL DPLINT(N,X,Y,C)
XX = r
CALL DPOLVL (NDER, XX, YFIT, YP, N, X, C, WORK, IERR)
intplt = YFIT
        
ENDFUNCTION

FUNCTION cintplt(dat,rs,r)
IMPLICIT NONE
DOUBLE COMPLEX,INTENT(IN)       ::dat(:)
DOUBLE PRECISION,INTENT(IN)     ::rs(:)
DOUBLE PRECISION,INTENT(IN)     ::r
DOUBLE COMPLEX                  ::cintplt
DOUBLE PRECISION                ::X(4),Y(4),C(4)
INTEGER,PARAMETER               ::N = 4
INTEGER                         ::NDER = 0,IERR
DOUBLE PRECISION                ::XX,YFIT,YP,WORK
INTEGER                         ::i

cintplt = dcmplx( &
        intplt(real(dat),rs,r), &
        intplt(imag(dat),rs,r))
        
ENDFUNCTION

ENDMODULE 
