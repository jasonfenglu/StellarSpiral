subroutine read_collective2d(q_loc,fnum)
use common_params
implicit none
integer::fnum
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,NVAR)
double precision::time
double precision,dimension(:,:),allocatable::density,momx,momy
#ifdef ADIABATIC
double precision,dimension(:,:),allocatable::ene
#endif
integer::ic1,ic2,ic3,ic4
character(len=8)::flnm
character(len=20)::dsetname

! write initial condition and basic info
ic1         = (fnum          -mod(fnum,1000))/1000
ic2         = (mod(fnum,1000)-mod(fnum, 100))/100
ic3         = (mod(fnum, 100)-mod(fnum,  10))/10
ic4         =  mod(fnum,  10)

flnm( 1:1 ) = 'M'
flnm( 2:2 ) = char(ichar('0')+ic1)
flnm( 3:3 ) = char(ichar('0')+ic2)
flnm( 4:4 ) = char(ichar('0')+ic3)
flnm( 5:5 ) = char(ichar('0')+ic4)
flnm( 6:8 ) = '.h5'

allocate(density(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))
allocate(momx(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))
allocate(momy(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))
#ifdef ADIABATIC
allocate(ene(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf))
#endif

dsetname='time'
call read1d(time,1,1,0,0,flnm,dsetname)
dsetname='density'
call read2d(density,ncell(1),ncell(2),ncell_loc(1),ncell_loc(2),ibuf,jbuf,2,2,flnm,dsetname)
dsetname='momx'
call read2d(momx,ncell(1),ncell(2),ncell_loc(1),ncell_loc(2),ibuf,jbuf,2,2,flnm,dsetname)
dsetname='momy'
call read2d(momy,ncell(1),ncell(2),ncell_loc(1),ncell_loc(2),ibuf,jbuf,2,2,flnm,dsetname)
#ifdef ADIABATIC
dsetname='ene'
call read2d(ene, ncell(1),ncell(2),ncell_loc(1),ncell_loc(2),ibuf,jbuf,2,2,flnm,dsetname)
#endif

t = time
q_loc(:,:,1) = density
q_loc(:,:,2) = momx
q_loc(:,:,3) = momy
#ifdef ADIABATIC
q_loc(:,:,4) = ene
#endif

deallocate(density)
deallocate(momx)
deallocate(momy)
#ifdef ADIABATIC
deallocate(ene)
#endif
end subroutine
