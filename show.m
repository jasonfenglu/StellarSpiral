clear all
close all

filename='density.h5'
den = hdf5read(filename,'density');
x = hdf5read(filename,'xcoord');
y = hdf5read(filename,'ycoord');
den = den';
imagesc(x,y,(den));
axis xy
axis equal
axis tight
colorbar
%caxis([-100 100])
