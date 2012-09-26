clear all
close all
filename='M0008.h5';
den = hdf5read(filename,'density');
x = hdf5read(filename,'x');
y = hdf5read(filename,'y');
imagesc(x,y,(den'));
axis xy
axis equal
axis tight
colorbar