%clear all
close all

i=0;
filename='../force.h5';
fx = hdf5read(filename,'stellarfx');
fy = hdf5read(filename,'stellarfy');
den = hdf5read(filename,'density');
x = hdf5read(filename,'x');
y = hdf5read(filename,'y');
fx = fx';
fy = fy';
f = (fx.^2+fy.^2).^0.5;
imagesc(x,y,(den));
title(sprintf('frame %d',i)) 
axis xy
axis equal
axis tight
colorbar
%caxis([-.03 .03])

pause(0.001);

