function last(index)

i =length(dir('../*.h5'))-2;

filename=sprintf('../M%04d.h5',i);
momx = hdf5read(filename,'momx');
momy = hdf5read(filename,'momy');
den  = hdf5read(filename,'density');
x = hdf5read(filename,'x');
y = hdf5read(filename,'y');

momx = momx';
den  = den' ;
vx = momx./den;
vy = momy./den;

switch (index)
    case 1
        imagesc(x,y,(den));
    case 2
        imagesc(x,y,vx);
    case 3
        imagesc(x,y,vy);
end

title(sprintf('frame %d',i))
axis xy
axis equal
axis tight
colorbar
%caxis([0 10])


