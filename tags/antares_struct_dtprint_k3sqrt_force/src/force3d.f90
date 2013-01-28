subroutine force3d(q_loc,fx,fy,fz)
use common_params
implicit none
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR)
double precision::fx(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3))
double precision::fy(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3))
double precision::fz(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3))
double precision,dimension(:,:,:),allocatable::den
! implement your external force here !!

fx = 0.d0
fy = 0.d0
fz = 0.d0

#if NDIM==3
#ifdef GRAVITY
 allocate(den(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3)))
 den= q_loc(1:ncell_loc(1),1:ncell_loc(2),1:ncell_loc(3),1)
 call selfgravity3d(den)
 fx = fx+sgx
 fy = fy+sgy
 fz = fz+sgz 
 deallocate(den)
#endif
#endif

end subroutine
