clear all
close all
loaddata

len=length(x);

fx=sgx(len/2+1-2:len-2-2,len/2+1-2);
fy=sgy(len/2+1-2:len-2-2,len/2+1-2);

[xx yy]=meshgrid(x(3:len-2),y(3:len-2));
den = transpose(den(3:len-2,3:len-2));
rr = sqrt(xx.^2+yy.^2);
dx = x(2)-x(1);
for k=1:len/2-2
  x0 = x(k+len/2);
  y0 = 0.5*(y(2)-y(1));
  rr0 = sqrt((xx-x0).^2+(yy-y0).^2);
  temp = GravConst*(xx-x0)./rr0.^3.*(den-153.0)*dx^2;
  ind = find(rr0 < 0.5*dx);
  temp(ind) = 0;
  ana_fx(k)=sum(temp(:));
end

for k=1:len/2-2
  x0 = x(k+len/2);
  y0 = 0.5*(y(2)-y(1));  
  rr0 = sqrt((xx-x0).^2+(yy-y0).^2);
  temp = GravConst*(yy-y0)./rr0.^3.*(den-153.0)*dx^2;
  ind = find(rr0 < 0.5*dx);
  temp(ind) = 0;
  ana_fy(k)=sum(temp(:));
end

figure(1)
plot(x(len/2+1:len-2),fx,'r+')
hold on
axis tight
xlabel('radius');
ylabel('x-force');
plot(x(len/2+1:len-2),ana_fx,'linewidth',2)
hold off

figure(2)
plot(y(len/2+1:len-2),fy,'ko')
hold on
axis tight
xlabel('radius');
ylabel('x-force');
plot(y(len/2+1:len-2),ana_fy,'linewidth',2)
hold off
