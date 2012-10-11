subroutine bndflux2d(q_loc,Fx,dd)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::Fx(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
integer::dd

! implement your boundary flux here~~~~~~


end subroutine
