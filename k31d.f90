program test
USE STELLARDISK
DOUBLE PRECISION        ::ri,rf,dr,r
DOUBLE PRECISION,ALLOCAtable ::dat(:,:),rs(:)
INTEGER                 ::n,i,j
type(typspiral)         ::spiral
DOUBLE COMPLEX          ::outk(6)
DOUBLE COMPLEX          ::t
!LOGICAL                 ::FLAG=.true.


ri = 0d0
rf = 15d0
N  = 2000

dr = (rf-ri)/dble(n)

ALLOCATE(dat(6,N))
ALLOCATE(rs(N))

CALL stdpara.readstd
CALL spiral.init(50,12.d0,stdpara,2)

!!start iteration
j = 0
DO WHILE(.true.)

IF(J<=14.and.J>=1)THEN
        write(*,*)'change to:'
        read(*,*)spiral.para(J)
ENDIF

CALL spiral.readw(2)
!CALL FindSpiral(spiral)
DO i = 1, N
       r = dr*dble(i-1)
       t = k3sqrt(r,spiral.ptr,outk)
       dat(1,i) = real(t)
       dat(2,i) = real(outk(1))
       dat(3:6,i) = real(outk(3:6))
       rs(i) = r
ENDDO

CALL plot
write(*,*)'choose which one to change.'
write(*,*)'      Lh        = this.para(1)'
write(*,*)'      rhoh      = this.para(2)'
write(*,*)'      Mb        = this.para(3)'
write(*,*)'      rb        = this.para(4)'
write(*,*)'      dM        = this.para(5)'
write(*,*)'      da        = this.para(6)'
write(*,*)'      db        = this.para(7)'
write(*,*)'      Qod       = this.para(8)'
write(*,*)'      q         = this.para(9)'
write(*,*)'      rq        = this.para(10)'
write(*,*)'      a1        = this.para(11)'
write(*,*)'      a2        = this.para(12)'
write(*,*)'      M1        = this.para(13)'
write(*,*)'      M2        = this.para(14)'
CALL spiral.print
read(*,*)j
if(j.eq.-1)EXIT


ENDDO
DEALLOCATE(dat)
DEALLOCATE(rs)
STOP
CONTAINS

SUBROUTINE plot
INTEGER                         ::PGOPEN
REAL                            ::TR(6)=0.
REAL                            ::zmax,zmin
REAL                            ::BRIGHT,CONTRA

IF(PGOPEN('/xs').NE.1)STOP
BRIGHT = 0.5
CONTRA = 0.9

CALL PGENV(0.,12.,-5.,2.,0,0)
DO i = 1, 6
        CALL PGSCI(i)
!       CALL PGLINE(N,real(rs),real(dat(i,:)),0)
        CALL PGPT(N,real(rs),real(dat(i,:)),0)
ENDDO

CALL PGCLOS
ENDSUBROUTINE

end      
