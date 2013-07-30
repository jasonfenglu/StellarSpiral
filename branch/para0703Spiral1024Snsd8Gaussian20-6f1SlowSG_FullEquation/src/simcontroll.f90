module simcontroll
IMPLICIT NONE
type simcon_type
        DOUBLE PRECISION                ::tend
        DOUBLE PRECISION                ::dtout
        INTEGER                         ::ncir
        INTEGER                         ::nframe
endtype
type(simcon_type),SAVE                  ::simcon
CONTAINS

SUBROUTINE      init_simcon(dtout,tend)
IMPLICIT NONE
DOUBLE PRECISION                        ::pi
DOUBLE PRECISION                        ::dtout,tend
DOUBLE PRECISION                        ::pspd
DOUBLE PRECISION                        ::W(4)
LOGICAL                                 ::bndu0(2)
namelist /SPIRALNML/w,bndu0
!$OMP CRITICAL
open(10,file='para.list')
read(10,nml=spiralnml)
close(10)
!$OMP END CRITICAL

pi = atan(1.d0)*4.d0

pspd = w(3)/2.d0
simcon%ncir   = 10
simcon%nframe = 800
simcon%tend   = 2.d0*pi/pspd*dble(simcon%ncir)
simcon%dtout  = simcon%tend/dble(simcon%nframe)

dtout = simcon%dtout
tend  = simcon%tend

ENDSUBROUTINE

endmodule
