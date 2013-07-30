clear all
close all

ffilename = '../force.h5';
spiral_den= hdf5read(ffilename,'density');
sx = hdf5read(ffilename,'x');
sy = hdf5read(ffilename,'y');
sx = sx(3:size(sx)-2);
sy = sy(3:size(sy)-2);

for i = 1:size(spiral_den,1)
    for j = 1:size(spiral_den,2)
        if(spiral_den(i,j)<0)
            spiral_den(i,j)=0;
        end
    end
end
spiral_den = spiral_den';


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

% [xx yy] = meshgrid(x,y);
% rr = sqrt(xx.^2+yy.^2);
% indr = find(rr > 7.3);
% den(indr) = NaN;
imagesc(x,y,(den));

% hold on;
% [C,h]=contour(sx,sy,spiral_den);
% ch = get(h,'child'); 
% alpha(ch,1);


title(sprintf('frame %d',i)) 
axis xy
axis equal
axis tight
colorbar
caxis([0 30])

pause(0.01);
end
