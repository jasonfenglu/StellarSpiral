PROGRAM RC
USE STELLARDISK
IMPLICIT NONE
INTEGER                         ::i
DOUBLE PRECISION                ::ri=0.d0
DOUBLE PRECISION                ::rf=20.d0
INTEGER                         ::N = 500
DOUBLE PRECISION                ::dr,r
DOUBLE PRECISION,ALLOCATABLE            ::dat(:,:)
DOUBLE PRECISION,ALLOCATABLE,TARGET     ::std(:,:)
DOUBLE PRECISION,TARGET         ::stdpara(10)
INTEGER                         ::PGBEG
!Halo
DOUBLE PRECISION                ::Lh,rhoh,gHalo,VHalo
!bulge
DOUBLE PRECISION                ::rb,Mb,gBulge,VBulge
!disk
DOUBLE PRECISION                ::dM,da,db,VDisk
!Toomre Q
DOUBLE PRECISION                ::Q,Qod,rq
INTEGER                         ::iter,j

!set up grid
dr = (rf-ri)/dble(N)
ALLOCATE(std(3,N))
do i = 1,N
        !r axis
        std(1,i) = ri+dr*dble(i)
enddo

Lh   = 2.8d0
rhoh = 4.0e7
Mb   = 10.0d7
rb   = 2.0d0
dM   = 7.0d10
da   = 2.7
db   = 0.3
Qod  = 1.d0
q    = 1.2d0
rq   = 2.8d0

stdpara = (/Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq/)
para=>stdpara
CALL INIT_STELLARDISK(N,rf)

!fill in standard data
do i = 1, N
        r        = std(1,i)
        std(2,i) = Omega(r)*r
        !std(3,i) = k3sqrt(r)
        do j = 1,4*N
                if(real(u(1,j)).gt.r)then
                        std(3,i) = &
                        real(h1(j))*sigma0(r)/snsd(r)**2/sigma0(r)
                        exit
                endif
        enddo
enddo
CALL ENDSTELLARDISK



!set pgplot
IF (PGBEG(0,'/xserve',1,1) .NE. 1) STOP
CALL PGSVP(0.3,0.70,0.3,0.70)
CALL PGSUBP(4,3)


!start interation to all components
ALLOCATE(dat(3,N))
do iter = 1,1 
        print *,'iter:',iter
        CALL iteration(iter,stdpara,dat,std,N)
enddo

CALL PGIDEN
CALL PGCLOS

DEALLOCATE(dat,std)
STOP
END PROGRAM

SUBROUTINE intitype(set_para)
USE STELLARDISK
DOUBLE PRECISION,TARGET         ::set_para(10)
para => set_para
ENDSUBROUTINE

SUBROUTINE plot(dat,std,N,num,iter)
use STELLARDISK,only:u,h1
IMPLICIT NONE
CHARACTER(len=10),DIMENSION(10)  ::lab
CHARACTER(len=20)               ::ch
DOUBLE PRECISION                ::dat(3,N),std(3,N)
DOUBLE PRECISION,PARAMETER      ::ri = 0.d0
DOUBLE PRECISION,PARAMETER      ::rf = 20.d0
DOUBLE PRECISION                ::num
INTEGER                         ::N,iter

lab(1) = 'Lh'
lab(2) = 'rhoh'
lab(3) = 'Mb ' 
lab(4) = 'rb ' 
lab(5) = 'dM ' 
lab(6) = 'da ' 
lab(7) = 'db ' 
lab(8) = 'Qod' 
lab(9) = 'q  ' 
lab(10)= 'rq ' 

!print *,'!!!!!!!!',lab(1)
write(ch,'(A6,E9.3E2)'),lab(iter),num
!print *,maxval(dat(3,:)),minval(dat(3,:))
CALL PGENV(REAL(ri),12.,-1e-2,1e-2,0,1)
CALL PGLINE(N,real(dat(1,:)),real(dat(3,:)))
CALL PGSLS(3)
CALL PGLINE(N,real(dat(1,:)),real(std(3,:)))
CALL PGSLS(1)
CALL PGLAB('','',ch)
!print *,std(3,:)-dat(3,:)
ENDSUBROUTINE

SUBROUTINE iteration(iter,stdpara,dat,std,N)
USE STELLARDISK
IMPLICIT NONE
!interface
!        SUBROUTINE iteration(iter,stdpara,dat,std,N)
!        IMPLICIT NONE
!        INTEGER                 ::iter
!        DOUBLE PRECISION,TARGET ::stdpara(10)
!        DOUBLE PRECISION,ALLOCATABLE            ::dat(:,:)
!        DOUBLE PRECISION,ALLOCATABLE,TARGET     ::std(:,:)
!        INTEGER                 ::N
!        ENDSUBROUTINE iteration
!endinterface
        
DOUBLE PRECISION                ::dat(3,N),std(3,N)
DOUBLE PRECISION                ::stdpara(10)
DOUBLE PRECISION,TARGET         ::sstdpara(10)
DOUBLE PRECISION,PARAMETER      ::div = 0.5d0
DOUBLE PRECISION                ::tmp,r
INTEGER                         ::iter,N,i,j,k,l

sstdpara = stdpara

j = 0
DO tmp = sstdpara(iter)*(1.d0-div),sstdpara(iter)*(1.d0+div),sstdpara(iter)*(div/5.d0)
!DO tmp = sstdpara(iter)*(1.d0),sstdpara(iter)*(1.d0)
!        print *,tmp
        sstdpara(iter) = tmp
        para=>sstdpara
        CALL INIT_STELLARDISK(N,20.d0)
        !fill in data
        do i = 1,N
                !r axis
                r = std(1,i)
                dat(1,i) = r
                dat(2,i) = Omega(r)*r
                !dat(3,i) = k3sqrt(r)
                
                do l = 1,4*N
                        if(real(u(1,l)).gt.r)then
                                dat(3,i) = &
                                real(h1(l))*sigma0(r)/snsd(r)**2/sigma0(r)
                                exit
                        endif
                enddo
        enddo
        !plot using pgplot
        CALL plot(dat,std,N,tmp,iter)
        print *,error()
        CALL ENDSTELLARDISK
        j = j + 1
!       print *,j
ENDDO
!padding to next page
DO k = j,11
        CALL PGPAGE
ENDDO

ENDSUBROUTINE
