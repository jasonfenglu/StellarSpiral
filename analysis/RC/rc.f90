PROGRAM RC
USE STELLARDISK
IMPLICIT NONE
INTEGER                         ::i
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::ri=0.d0
DOUBLE PRECISION                ::rf=20.d0
INTEGER                         ::N = 1000
DOUBLE PRECISION                ::dr,r
DOUBLE PRECISION,ALLOCATABLE    ::dat(:,:)
INTEGER                         ::PGBEG
CHARACTER                       ::cmd


!CALL getarg(1,arg)
!READ(arg,*)wr
!CALL getarg(2,arg)
!READ(arg,*)wi


!set up grid
dr = (rf-ri)/dble(N)
ALLOCATE(dat(5,N))
do i = 1,N
        !r axis
        dat(1,i) = ri+dr*dble(i)
enddo

!set basic parameters
QParameter=QParameterType(1.d0,3.9d0,3.1d0)
wr        = 26.d0
wi        = 0.d0

do while(.true.)
!fill in data
do i = 1, N
        !fill Omega
        r        = dat(1,i)
        dat(2,i) = Omega(r)
        dat(3,i) = Q(r)
        dat(4,i) = k3sqrt(r)
        dat(5,i) = nu(r)
enddo

!start pgplot
!IF (PGBEG(0,'?',1,1) .NE. 1) STOP
IF (PGBEG(0,'/xserve',1,1) .NE. 1) STOP
CALL PGSVP(0.0,0.95,0.0,0.95)
CALL PGSUBP(2,2)
CALL PGENV(REAL(ri),REAL(rf),0.,real(maxval(dat(2,:))),0,1)
!plot Omega
CALL PGLINE(N,real(dat(1,:)),real(dat(2,:)))
CALL PGLAB("","","Omega")
CALL PGENV(REAL(ri),REAL(rf),0.,real(maxval(dat(3,:))),0,1)
CALL PGLINE(N,real(dat(1,:)),real(dat(3,:)))
CALL PGLAB("","","Q")
CALL PGENV(REAL(ri),REAL(rf),-1.,5.,0,1)
CALL PGLINE(N,real(dat(1,:)),real(dat(4,:)))
CALL PGLAB("","","k3sqr")
CALL PGENV(REAL(ri),REAL(rf),-2.,2.,0,1)
CALL PGLINE(N,real(dat(1,:)),real(dat(5,:)))
CALL PGLAB("","","k3sqr")
CALL PGCLOS
write(*,*)QParameter
read(*,*)QParameter%Qod,QParameter%q,QParameter%rq
enddo


DEALLOCATE(dat)
STOP
END PROGRAM

