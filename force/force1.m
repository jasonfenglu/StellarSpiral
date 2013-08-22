clear all
close all

ibuf = 2;

filename='force.h5';
fx = hdf5read(filename,'fx');
fy = hdf5read(filename,'fy');
x = hdf5read(filename,'xcoord');
y = hdf5read(filename,'ycoord');
fa  = hdf5read(filename,'fa');

x = x(ibuf+1:length(x)-ibuf);
y = y(ibuf+1:length(y)-ibuf);


fx = fx';
fy = fy';

% fx = fx/(4.3e-6);
% fy = fy/(4.3e-6);

ff  = (fx.^2+fy.^2).^0.5;


r = x(length(x)/2+1:length(x));
f = ff(length(x)/2,length(x)/2+1:length(x));
faa = fa(length(x)/2,length(x)/2+1:length(x));
plot(r,f);
hold on;
plot(r,2*r.*exp(-r.^2),'r');
plot(r,faa,'o');

filename='antforce.h5';
fx = hdf5read(filename,'stellarfx');
fy = hdf5read(filename,'stellarfy');
ffa  = (fx.^2+fy.^2).^0.5;
r = x(length(x)/2+1:length(x));
f = ffa(length(x)/2,length(x)/2+1:length(x));
plot(r,f,'+');

figure;
imagesc(x,y,ff-ffa);
axis equal;
axis tight;
colorbar
%% 2d to 1d, not using
% r=[];
% f=[];
% 
% for i = 1:3:length(x)
%     for j = 1:3:length(y)
% 
%         fxx = fx(i,j);
%         fyy = fy(i,j);
%         rr  = (x(i)^2+y(j)^2)^.5;
%         ff  = (fxx^2+fyy^2)^.5;
%         r   = [r rr];
%         f   = [f ff];
%     end
% end
% 
% [r r1]=sort(r);
% ff = zeros(length(f));
% for i=1:length(f)
%     ff(i) = f(r1(i));
% end
% plot(r,f);
% hold on;
% plot(r,2*r.*exp(-r.^2));