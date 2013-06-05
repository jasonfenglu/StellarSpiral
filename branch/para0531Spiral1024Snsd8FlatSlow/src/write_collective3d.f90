subroutine write_collective3d(q_loc,fnum)
use common_params
implicit none
integer::fnum
double precision::q_loc(1-ibuf:ncell_loc(1)+ibuf,1-jbuf:ncell_loc(2)+jbuf,1-kbuf:ncell_loc(3)+kbuf,NVAR)
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

dsetname='x'
call output1d(x_loc,ncell(1),ncell_loc(1),ibuf,2,flnm,dsetname,1) !1: new 0:open
dsetname='y'
call output1d(y_loc,ncell(2),ncell_loc(2),jbuf,2,flnm,dsetname,0)
dsetname='z'
call output1d(z_loc,ncell(3),ncell_loc(3),kbuf,2,flnm,dsetname,0)
dsetname='time'
call output1d(t,1,1,0,0,flnm,dsetname,0)
dsetname='xrange'
call output1d(xrange,2,2,0,0,flnm,dsetname,0)
dsetname='yrange'
call output1d(yrange,2,2,0,0,flnm,dsetname,0)
dsetname='ncell'
call output1d(dble(ncell),3,3,0,0,flnm,dsetname,0)
dsetname='snd'
call output1d(snd,1,1,0,0,flnm,dsetname,0)
dsetname='density'
call output3d(q_loc(:,:,:,1),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
dsetname='momx'
call output3d(q_loc(:,:,:,2),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
dsetname='momy'
call output3d(q_loc(:,:,:,3),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
dsetname='momz'
call output3d(q_loc(:,:,:,4),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)

#ifdef ADIABATIC
#ifdef MHD
dsetname='bxl'
call output3d(q_loc(:,:,:,5),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
dsetname='byl'
call output3d(q_loc(:,:,:,6),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
dsetname='bzl'
call output3d(q_loc(:,:,:,7),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
dsetname='ene'
call output3d(q_loc(:,:,:,8),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
dsetname='bxr'
call output3d(q_loc(:,:,:,9),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
dsetname='byr'
call output3d(q_loc(:,:,:,10),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
dsetname='bzr'
call output3d(q_loc(:,:,:,11),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)

#else
dsetname='ene'
call output3d(q_loc(:,:,:,5),ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),ibuf,jbuf,kbuf,2,2,2,flnm,dsetname,0)
#endif
#endif




end subroutine
