close all
clear all
flnm='M0001.h5'
x=hdf5read(flnm,'x');
y=hdf5read(flnm,'y');
t=hdf5read(flnm,'time');
ncell=hdf5read(flnm,'ncell');
xrange=hdf5read(flnm,'xrange');
yrange=hdf5read(flnm,'yrange');
snd=hdf5read(flnm,'snd');
den=(hdf5read(flnm,'density'));
momx=(hdf5read(flnm,'momx'));
momy=(hdf5read(flnm,'momy'));
sgx=(hdf5read(flnm,'sgx'));
sgy=(hdf5read(flnm,'sgy'));
%den=transpose(hdf5read(flnm,'density'));
%momx=transpose(hdf5read(flnm,'momx'));
%momy=transpose(hdf5read(flnm,'momy'));
%sgx=transpose(hdf5read(flnm,'sgx'));
%sgy=transpose(hdf5read(flnm,'sgy'));
GravConst=hdf5read(flnm,'GravConst');
