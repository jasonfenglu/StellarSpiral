clear all
close all

for i = 0:length(dir('../*.h5'))-1
filename=sprintf('../M%04d.h5',i);
den = hdf5read(filename,'density');
px  = hdf5read(filename,'momx');
py  = hdf5read(filename,'momy');
x = hdf5read(filename,'x');
y = hdf5read(filename,'y');
den = den';
px  = px';
py  = py';

vx  = px./den;
vy  = py./den;

imagesc(x,y,(den));
title(sprintf('frame %d',i)) 
axis xy
axis equal
axis tight
colorbar
%caxis([-1000 1000])

pause(0.01);
end
