%clear all
close all

filename='../force.h5';
fx = hdf5read(filename,'stellarfx');
fy = hdf5read(filename,'stellarfy');
den = hdf5read(filename,'density');
x = hdf5read(filename,'x');
y = hdf5read(filename,'y');
fx = fx';
fy = fy';

x = x(3:length(x)-4);
y = y(3:length(y)-4);

r=[];
f=[];

for i = 1:3:length(x)
    for j = 1:3:length(y)

        fxx = fx(i,j);
        fyy = fy(i,j);
        rr  = (x(i)^2+y(j)^2)^.5;
        ff  = (fxx^2+fyy^2)^.5;
        r   = [r rr];
        f   = [f ff];
    end
end

[r r1]=sort(r);
ff = zeros(length(f));
for i=1:length(f)
    ff(i) = f(r1(i));
end
plot(r,f);
axis([6.5 7 0 0.03]);