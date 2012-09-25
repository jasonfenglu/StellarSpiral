subroutine force2dCTU(q_loc,fx,fy)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::fx(1:ncell_loc(1),1:ncell_loc(2)),fy(1:ncell_loc(1),1:ncell_loc(2))
double precision,dimension(:,:),allocatable::den
! implement your external force here !!

fx = 0.d0
fy = 0.d0

#if NDIM==2
#ifdef GRAVITY
 allocate(den(1:ncell_loc(1),1:ncell_loc(2)))
 den= q_loc(1:ncell_loc(1),1:ncell_loc(2),1)
 call selfgravity2d(den)
 fx = fx+sgx
 fy = fy+sgy 
 deallocate(den)
#endif
#endif

end subroutine
