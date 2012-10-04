        MODULE density
        IMPLICIT NONE
        REAL*8                  ::alp
        REAL*8                  ::alpcsq2
        CHARACTER(len=4)        ::density_type
        REAL*8                  ::density_max
        SAVE                    alp,alpcsq2,density_type
        CONTAINS

        SUBROUTINE init_density(option,sspd,dmax)
        IMPLICIT NONE
        CHARACTER(*)            ::option
        REAL*8                  ::sspd
        REAL*8                  ::dmax
        if(option .eq. 'bell')then
              alp = 1.d0/48.8581d0
              alpcsq2 = 2.d0 * alp * sspd **2.d0

        elseif(option .eq. 'flat')then
              alp = 0.d0
              alpcsq2 = 0.d0

        endif
        density_type = option
        density_max  = dmax
        ENDSUBROUTINE

        FUNCTION density_distribution(r)
        IMPLICIT NONE
        REAL*8                  ::r,max_density
        REAL*8                  ::density_distribution

        !write(*,*)alp
        
        if(density_type .eq. 'bell')then
                density_distribution = 
     &          density_max*dexp(-alp*r**2)
        elseif(density_type .eq. 'flat')then
                density_distribution = density_max
        endif
        ENDFUNCTION density_distribution
        

        ENDMODULE density
