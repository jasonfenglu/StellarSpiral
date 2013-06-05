clear all
close all

filename='force.h5';
f = hdf5read(filename,'f');
x = hdf5read(filename,'xcoord');
y = hdf5read(filename,'ycoord');
r = hdf5read(filename,'r');
fr = hdf5read(filename,'fr');
fc = hdf5read(filename,'centri');

f = f';

imagesc(x,y,(f));
axis xy
axis equal
axis tight
colorbar

figure
plot(r,fr)

%caxis([-100 100])

% 
% filename='antforce.h5'
% fx = hdf5read(filename,'stellarfx');
% fy = hdf5read(filename,'stellarfy');
% f  = (fx.^2+fy.^2).^.5;
% f = f';
% x2 = hdf5read(filename,'x');
% y2 = hdf5read(filename,'y');
% x2 = x2(3:length(x2)-2);
% y2 = y2(3:length(y2)-2);
% 
% imagesc(x,y,den-f);
% axis xy
% axis equal
% axis tight
% colorbar