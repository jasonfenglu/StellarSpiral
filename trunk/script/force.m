clear all
close all

for i = 0:length(dir('../*.h5'))-1
filename=sprintf('../M%04d.h5',i);
fx = hdf5read(filename,'forcex');
fy = hdf5read(filename,'forcey');
x = hdf5read(filename,'x');
y = hdf5read(filename,'y');
fx = fx';
fy = fy';
%fx = (fx.^2+fy.^2).^0.5;
imagesc(x,y,(fx));
title(sprintf('frame %d',i)) 
axis xy
axis equal
axis tight
colorbar
%caxis([-80 80])

pause(0.001);
end
