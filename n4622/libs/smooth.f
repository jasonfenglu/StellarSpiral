!$Id: smooth.f 118 2011-07-21 05:07:37Z ccfeng $
        MODULE smooth
        IMPLICIT NONE
        REAL*8                  ::pi = datan(1.d0)*4.d0

        interface
                !FUNCTION exp_smooth(r,r0)
                !REAL*8          ::r
                !REAL*8,OPTIONAL ::r0
                !ENDFUNCTION exp_smooth
        END interface

        CONTAINS

        FUNCTION cosine_smooth(r)
        IMPLICIT NONE
        REAL*8                  ::cosine_smooth
        REAL*8                  ::r
        REAL*8                  ::r0 = 15.d0

        if (r>r0) then
        cosine_smooth = 
     &  5.d-1*(dcos((r-r0)*pi/5.d0)+1.d0)
        
        else
        cosine_smooth = 1.d0
        ENDIF

        ENDFUNCTION cosine_smooth

        FUNCTION exp_smooth(r,start)
        REAL*8                  ::exp_smooth
        REAL*8                  ::r
        REAL*8,optional         ::start 

        REAL*8                  ::r0

        if(present(start))then
                r0 = start
        else
                r0 = 15.d0
        endif

        if (r>r0) then
        exp_smooth = 
     &   dexp(-(r-r0)**2/12.d0)

        else
        exp_smooth =  1.d0
        ENDIF

        ENDFUNCTION exp_smooth

        FUNCTION linear_increase(t,tf)
        IMPLICIT NONE
        REAL*8                  ::linear_increase
        REAL*8                  ::t,tf
        
        linear_increase = dmin1(t/tf*10.d0,1.d0)
        ENDFUNCTION linear_increase

        ENDMODULE smooth
