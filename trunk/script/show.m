clear all
close all

for i = 0:length(dir('../*.h5'))-1
filename=sprintf('../M%04d.h5',i);
den = hdf5read(filename,'density');
x = hdf5read(filename,'x');
y = hdf5read(filename,'y');
den = den';
imagesc(x,y,(den));
title(sprintf('frame %d',i)) 
axis xy
axis equal
axis tight
colorbar
%caxis([-100 100])

pause(0.001);
end
