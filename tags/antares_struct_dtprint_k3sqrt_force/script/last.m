clear all
close all

i =length(dir('../*.h5'))-2;
filename=sprintf('../M%04d.h5',i);
momx = hdf5read(filename,'momx');
den  = hdf5read(filename,'density');
x = hdf5read(filename,'x');
y = hdf5read(filename,'y');
den = den'./momx';
imagesc(x,y,(den));
title(sprintf('frame %d',i))
axis xy
axis equal
axis tight
colorbar
%caxis([0 10])


