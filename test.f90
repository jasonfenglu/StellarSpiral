program test
USE STELLARDISK
DOUBLE PRECISION        ::ri,rf,dr,r
DOUBLE COMPLEX          ::h1,t
INTEGER                 ::n,i
type(typspiral)spiral
type tst
     integer,allocatable::dat(:)
endtype
type(tst) tt
DOUBLE COMPLEX          ::Y(2,100),Y1(2)
DOUBLE PRECISION        ::tmp(2,2)
CHARACTER(len=225)      ::CH
DOUBLE PRECISION        ::a,b
LOGICAL                 ::FLAG
DOUBLE COMPLEX          ::outk(6)

!CALL GETCWD(CH)
!CH = TRIM(CH)
!print *,INDEX(CH,'/',.true.)
!CH = TRIM(CH(INDEX(CH,'/',.true.)+1:LEN(CH)))
!print *,CH,LEN(CH)
!print *,TRIM(CH)
!STOP

ri = 0.d0
rf = 1.3d1
!rf = 10.d0
N  = 100

dr = (rf-ri)/dble(n)

!tmp = 2.d0
!CALL h5write(tmp,2,2,"dsetf.h5","dset")
!CALL h5write(tmp(1,:),2,"dsetf.h5","dset2")
!CALL h5read(tmp,2,2,"dsetf.h5","dset")
!print *,tmp

!!test file names stack
!write(6,*)'return',checkfileext('123')
!write(6,*)'return',checkfileext('234')
!write(6,*)'return',checkfileext('123')
!call printpt

!CALL stdpara.printpara

CALL stdpara.readstd
CALL spiral.init(100,15.d0,stdpara,2)
CALL spiral.readw(2)
CALL FindSpiral(spiral.ptr)
!CALL spiral.free
!print *,real(spiral.u)
!DO i = 1, 100
!        r = dr*dble(i-1)
!       !write(6,*)r,kappa(r,Spiral)*8.d0/pi/g/gasdensity(r,0.d0)
!       !write(6,*)r,abs(k3sqrt(r,spiral)),spiral.co
!       t = k3sqrt(r,spiral.ptr,outk)
!       !print *,r,snsd(r,spiral),ToomreQ(r,spiral)
!       !write(*,'(4(D15.5))')r,real(t),imag(t),abs(t)
!!      print *,real(spiral.u(1,i)),real(spiral.h1(i)),imag(spiral.h1(i))
!        write(6,'(7(D15.5))')r,kappa(r,spiral.ptr)
!ENDDO


!CALL stdpara.readstd
!ALLOCATE(tt.dat(3))
!!$OMP PARALLEL FIRSTPRIVATE(spiral)
!!$OMP DO 
!DO i = 1, 32
!CALL spiral.init(spiral,500,12.d0,stdpara,1)
!CALL FindSpiral(spiral)
!CALL spiral.final
!print *,spiral.error
!ENDDO
!!$OMP END PARALLEL

!CALL stdpara.readstd
!CALL spiral.init(N,12.d0,stdpara,1)
!CALL spiral.readw(1)
!CALL FindSpiral(spiral)
!!CALL FindPhi1(spiral)
!DO i = 1, N
!        r = dr*dble(i)
!!      write(6,*)r,real(spiral.u(2,i)),imag(spiral.u(2,i))
!!      write(6,*)r,real(k3sqrt(r,spiral)),imag(k3sqrt(r,spiral))
!!       write(6,*)r,real(cintplt(spiral.h1,spiral.r,r)),&
!!       real(cintplt(spiral.phi1r,spiral.r,r))
!        write(6,*)r,abs(cintplt(spiral.h1,spiral.r,r))/snsd(r,spiral)**2,sigma0(r,spiral)
!ENDDO

!CALL stdpara.readstd
!CALL spiral.init(500,12.d0,stdpara,2)
!CALL spiral.readw(2)
!CALL FindSpiral(spiral)
!DO i = 1, N
!        r = dr*dble(i)
!!       print *,spiral.r(i),&
!!       real(spiral.h1(i)),intplt(real(spiral.h1),spiral.r,r)
!!       print *,r,intplt(real(spiral.h1),spiral.r,r)
!        print *,r,sigma1r(r,spiral)
!ENDDO

!Y1(1) = dcmplx(0.d0,0.d0)
!Y1(2) = dcmplx(1.d0,0.d0)
!CALL rk4d1(0.d0,3.d0,100,F,F3,Y,Y1)
!DO i = 1, 100
!        print *,real(Y(1,i)),Y(2,i)
!ENDDO
!
!stop
CONTAINS


!FUNCTION F(r)
!IMPLICIT NONE
!DOUBLE PRECISION                ::r
!DOUBLE PRECISION                ::F
!F = sin(r)
!ENDFUNCTION
!
!Function F2(r) result(ans)
!USE NUM
!        IMPLICIT NONE
!        DOUBLE PRECISION,INTENT(IN)     ::r
!        DOUBLE PRECISION                ::a,rr
!        DOUBLE PRECISION                ::ans
!        DOUBLE PRECISION                ::ferr = 1.d-15
!        INTEGER                         ::IERR,K=6000
!        a = zerolimit
!        rr = r
!        IERR = 0
!        if(rr.eq.0.d0)then
!                ans = 0.d0
!        else
!                CALL DGAUS8(F,a,rr,fERR,ans,IERR)
!!               CALL DQNC79(F,a,rr,fERR,ans,IERR,K)
!        endif
!
!ENDFUNCTION 
!
!FUNCTION F3(r)
!IMPLICIT NONE
!DOUBLE PRECISION                ::r
!DOUBLE COMPLEX                  ::F3
!F3 = (1.d0,1.d0)*r**2
!ENDFUNCTION

end      
