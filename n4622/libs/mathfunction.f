!$Id: mathfunction.f 158 2011-08-10 06:57:09Z ccfeng $       
        MODULE MATH

        CONTAINS
        FUNCTION ellipticK(k)
        !DEC$ ATTRIBUTES ALIAS : 'wrapper_ellipticK'::wrapper_ellipticK
        IMPLICIT NONE
        REAL*8          ::k
        REAL*8          ::ellipticK
        REAL*8          ::ans

        call wrapper_ellipticK(k,ans)
        ellipticK = ans
        ENDFUNCTION ellipticK

        FUNCTION ellipticE(k)
        !DEC$ ATTRIBUTES ALIAS : 'wrapper_ellipticE'::wrapper_ellipticE
        IMPLICIT NONE
        REAL*8          ::k
        REAL*8          ::ellipticE
        REAL*8          ::ans

        call wrapper_ellipticE(k,ans)
        ellipticE = ans
        ENDFUNCTION ellipticE
        END MODULE MATH
