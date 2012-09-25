module common_params

integer::ibuf=2,jbuf=2,kbuf=2 !number of ghost zone
double precision::snd  ! sound speed
double precision:: xrange(2),yrange(2),zrange(2)
double precision::tend,CFL,dtout ! final stopping time, Courant number, time step for data output
integer::ncell(3),ncell_loc(3),nprocs,myid ! global cell number, local cell number, number of processor
double precision,dimension(:),allocatable::x_loc,y_loc,z_loc
double precision,dimension(:),allocatable::x,y,z
double precision:: dx,dy,dz,dt,t ! cell size in x,y,z, global time step, global time
double precision:: gam,GravConst
#if NDIM==2
#ifdef CTU
integer,parameter:: NVAR=3
#else
#ifdef MHD
integer,parameter:: NVAR=11
#else
integer,parameter:: NVAR=4
#endif
#endif
#endif

#if NDIM==3
#ifdef MHD
integer,parameter:: NVAR=11
#else
integer,parameter:: NVAR=5
#endif
#endif

#if NDIM==2
double precision,dimension(:,:),allocatable::sgx,sgy !selfgravity
double complex  ,dimension(:,:),allocatable::sgxker,sgyker !kernal 
#endif
#if NDIM==3
double precision,dimension(:,:,:),allocatable::sgx,sgy,sgz 
double complex  ,dimension(:,:,:),allocatable::sgxker,sgyker,sgzker
#endif

end module common_params
