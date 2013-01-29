PROGRAM RC
USE STELLARDISK
IMPLICIT NONE
INTEGER                         ::i
CHARACTER(len=32)               ::arg
DOUBLE PRECISION                ::ri=0.d0
DOUBLE PRECISION                ::rf=20.d0
INTEGER                         ::N = 200
DOUBLE PRECISION                ::dr,r
DOUBLE PRECISION,ALLOCATABLE    ::dat(:,:)
REAL                            ::LoweRC(20,2),tmp
INTEGER                         ::PGBEG
CHARACTER                       ::cmd


!CALL getarg(1,arg)
!READ(arg,*)wr
!CALL getarg(2,arg)
!READ(arg,*)wi

open(unit=10,file='RC_Lowe')
do i = 1,20
        read(10,*)LoweRC(i,1),LoweRC(i,2)
enddo

!set up grid
dr = (rf-ri)/dble(N)
ALLOCATE(dat(8,N))
do i = 1,N
        !r axis
        dat(1,i) = ri+dr*dble(i)
enddo

!set basic parameters
QParameter=QParameterType(1.d0,3.9d0,3.1d0)
wr        = 26.d0
wi        = 0.d0

!fill in data
do i = 1, N
        !fill Omega
        r        = dat(1,i)
        dat(2,i) = Omega(r)*r
        dat(5,i) = Omega(r)-kappa(r)/2.d0
        dat(6,i) = Sigma0(r)
enddo

do while(.true.)
do i = 1, N
        dat(3,i) = Q(r)
        dat(4,i) = k3sqrt(r)

enddo
!start pgplot
!IF (PGBEG(0,'?',1,1) .NE. 1) STOP
IF (PGBEG(0,'/xserve',1,1) .NE. 1) STOP
CALL PGSVP(0.0,0.95,0.0,0.95)
CALL PGSUBP(2,2)
!plot Velocity
tmp  = amax1(maxval(LoweRC(:,2)),real(maxval(dat(2,:))))*1.1
CALL PGENV(REAL(ri),REAL(rf),0.,tmp,0,1)
CALL PGLINE(N,real(dat(1,:)),real(dat(2,:)))
CALL PGSCI(2)
CALL PGPT(20,LoweRC(:,1),LoweRC(:,2),0)
CALL PGSCI(1)
CALL PGLAB("","","V")
!plot Q
CALL PGENV(REAL(ri),REAL(rf),0.,real(maxval(dat(3,:))),0,1)
CALL PGLINE(N,real(dat(1,:)),real(dat(3,:)))
CALL PGLAB("","","Q")
!plot k3sqr
CALL PGENV(REAL(ri),REAL(rf),-0.1,1.,0,1)
CALL PGLINE(N,real(dat(1,:)),real(dat(4,:)))
CALL PGLAB("","","k3sqr")
!plot nu
CALL PGENV(REAL(ri),REAL(rf),0.,real(maxval(dat(6,:))),0,1)
CALL PGLINE(N,real(dat(1,:)),real(dat(6,:)))
CALL PGLAB("","","snsd")
CALL PGIDEN
CALL PGCLOS
write(*,*)QParameter
read(*,*)QParameter%Qod,QParameter%q,QParameter%rq
enddo


DEALLOCATE(dat)
STOP
END PROGRAM

