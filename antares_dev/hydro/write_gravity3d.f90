subroutine write_gravity3d(fnum)
use common_params
implicit none
integer::fnum
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

dsetname='GravConst'
call output1d(GravConst,1,1,0,0,flnm,dsetname,0)
dsetname='sgx'
call output3d(sgx,ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),0,0,0,0,0,0,flnm,dsetname,0)
dsetname='sgy'
call output3d(sgy,ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),0,0,0,0,0,0,flnm,dsetname,0)
dsetname='sgz'
call output3d(sgz,ncell(1),ncell(2),ncell(3),ncell_loc(1),ncell_loc(2),ncell_loc(3),0,0,0,0,0,0,flnm,dsetname,0)

end subroutine
