        MODULE slatec
        real*8                :: beta
        SAVE beta
        CONTAINS

        FUNCTION laplace_coefi(m,beta)
        implicit none
        real*8                :: laplace_coefi,m,beta
        real*8                :: h,s                 ! the form of the formula
        h = 0.d0
        s = 0.5d0
        laplace_coefi = dqawoe_laplace_coefi(m,beta,h,s)
        ENDFUNCTION laplace_coefi

        FUNCTION dlaplace_coefi(m,beta)
        implicit none
        real*8                :: dlaplace_coefi,m,beta
        real*8                :: h,s                 ! the form of the formula

        h = 1.d0
        s = 1.5d0
        dlaplace_coefi = dqawoe_laplace_coefi(m,beta,h,s)
        ENDFUNCTION dlaplace_coefi

        FUNCTION dqawoe_laplace_coefi(m,in_beta,h,s)
        implicit none
        real*8                :: dqawoe_laplace_coefi,in_m,in_beta
        !real*8,external       :: f
        real*8                :: a = 0.d0
        real*8                :: b = datan(1.d0)*4.d0
        real*8                :: omega
        real*8                :: m
        real*8                :: h,s                 ! the form of the formula
        integer               :: integr = 1
        real*8                :: epsabs = 10.d-12
        real*8                :: epsrel = 10.d-12
        integer               :: limit = 1000
        integer               :: icall = 1
        integer               :: maxp1 = 5
        real*8                :: result
        real*8                :: abserr
        integer               :: neval
        integer               :: ier
        integer               :: last
        real*8,ALLOCATABLE    :: alist(:),blist(:)
        real*8,ALLOCATABLE    :: rlist(:),elist(:)
        integer,ALLOCATABLE   :: iord(:), nnlog(:)
        integer               :: momcom 
        real*8,ALLOCATABLE    :: chebmo(:,:)


        real*8 pi

        beta = in_beta

        !Allocate:
        ALLOCATE(alist(limit))
        ALLOCATE(blist(limit))
        ALLOCATE(rlist(limit))
        ALLOCATE(elist(limit))
        ALLOCATE(iord(limit))
        ALLOCATE(nnlog(limit))
        ALLOCATE(chebmo(maxp1,25))

        pi = datan(1.d0)*4.d0
        omega = m

        call dqawoe (f,a,b,omega,integr,epsabs,epsrel,limit,icall,
     *               maxp1,result,abserr,neval,ier,last,alist,blist,
     *               rlist,elist,iord,nnlog,momcom,chebmo)

c       write( *, '(a,g14.6)' )'Integral left endpoint A =    ', a
c       write( *, '(a,g14.6)' )'Integral right endpoint B =   ', b
c       write( *, '(a,g14.6)' )'Estimated integral is         ',result
c       write( *, '(a,g14.6)' )'Estimated integral error =    ',abserr
        
        if (IER .GT. 0)then
                write(*,*)'dqawoe integral fail'
                stop
        endif
        dqawoe_laplace_coefi = result

        contains
        function f(x)
        implicit none
        real*8                ::f
        real*8                ::x
        real*8                ::pi
        pi = datan(1.d0)*4.d0


        f = 2.0/PI*((dcos(x)-beta)**h)
     c    /(1+beta**2-2*beta*dcos(x))**s

        return
        endfunction f

        ENDFUNCTION dqawoe_laplace_coefi

        ENDMODULE slatec
