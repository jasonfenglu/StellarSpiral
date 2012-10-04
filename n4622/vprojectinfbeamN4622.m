%$Id: vprojectinfbeamN4622.m 112 2011-07-20 03:18:55Z ccfeng $
% clear all;
Nx = 516;
%close all;
load M0500;
q       = reshape(q,Nx,Nx,3);
rho     = q(:,:,1);
ux  = q(:,:,2)./rho;
uy  = q(:,:,3)./rho;
% Nx  = length(rho);
Ny  = Nx;
range = 30;
buffer = 2;
dx=range/(Nx-2*buffer);
dy=dx;

phi = (10+90)/180*pi;  % position angle
alpha =(112)/180*pi;  % position angle of line of node
theta =-19/180*pi;  % inclination angle

center = length(rho)/2;
a = -center:1:length(rho)-center-1;
b = -center:1:length(rho)-center-1;
[x,y] = meshgrid(a,b);
xx = x*cos(phi)-y*sin(phi);
yy = x*sin(phi)+y*cos(phi);
uxx = ux*cos(phi)-uy*sin(phi);
uyy = ux*sin(phi)+uy*cos(phi);

xxx = cos(alpha)*(xx*cos(alpha)+yy*sin(alpha))*(1-cos(theta))+xx*cos(theta);
yyy = sin(alpha)*(xx*cos(alpha)+yy*sin(alpha))*(1-cos(theta))+yy*cos(theta);
zzz = yy*cos(alpha)*sin(theta)-xx*sin(alpha)*sin(theta);
 
uzzz = uyy*cos(alpha)*sin(theta)-uxx*sin(alpha)*sin(theta);


%%%%% beam begin %%%%
%sigma=0.0768*4.5;%velocity field
sigma=0.0768*4.5;%velocity field
% sigma2=0.0504/3;%velocity field
% ang = 2.77/180*pi;%velocity field
% sigma1=0.02679/3;%1.3mm continuum
% sigma2=0.05833/3;%1.3mm continuum
% ang = -0.8/180*pi;%1.3mm continuum


exuzzz=zeros(2*Nx);
exuzzz(1:Nx, 1:Nx)=uzzz(1:Nx, 1:Nx);
aa=dx*(a+0.5);
bb=dy*(b+0.5);
[aaa,bbb]=meshgrid(aa,bb);
% raaa = aaa*cos(ang) + bbb*sin(ang);
% rbbb = -aaa*sin(ang) + bbb*cos(ang);
green=exp(-(aaa.^2+bbb.^2)/(2*sigma^2))/(2*pi*sigma*sigma);
figure
mesh(aaa,bbb,(green),'FaceColor','interp');
axis([-15 15 -15 15]);
axis square;
exgreen=zeros(2*Nx);
exgreen(1:Nx, 1:Nx)=green(1:Nx, 1:Nx);
Fexuzzz=fft2(exuzzz);
Fexgreen=fft2(exgreen);
Fresult=Fexuzzz.*Fexgreen;
result=ifft2(Fresult);
%result = conv2(uzzz, green)
exrho=zeros(2*Nx);
exrho(1:Nx, 1:Nx)=rho(1:Nx, 1:Nx);
Fexrho=fft2(exrho);
Fresultrho = Fexrho.*Fexgreen;
resultrho = ifft2(Fresultrho);
%mesh(xxx,yyy,transpose(result(514:1541, 514:1541)),'FaceColor','interp');
 xxxx = xxx*dx;
 yyyy = yyy*dy;
result2 = result(258:773, 258:773);
figure;
contour(yyyy,xxxx,(result2*dx*dy),20,'LineWidth',1);
axis square;
print('-depsc','convelocity.eps');
% %contour(xxxx,yyyy,(result2*dx*dy),20,'LineWidth',1,'Color','b')
% %contourf(xxxx,yyyy,(result2*dx*dy),20,'LineWidth',1)
figure;
contour(yyyy,xxxx,uzzz,20);
axis([-15 15 -15 15]);
axis square;
print('-depsc','velocity.eps');
% axis([-0.3 0.3 -0.3 0.3]);
% axis square;
% set(gca,'XTick',[-0.2 -0.1 0 0.1 0.2]);
% set(gca,'YTick',[-0.2 -0.1 0 0.1 0.2]);
% colorbar;
% %hold on;
% 
% resultrho2 = resultrho(258:773,258:773);
% figure(3)
% mesh(xxxx,yyyy,log10(resultrho2*dx*dy),'FaceColor','interp');
% %mesh(xxxx,yyyy,(log10(rho)),'FaceColor','interp');
% caxis([1.3 3]);
% axis([-0.3 0.3 -0.3 0.3]);
% %axis([-0.45 0.45 -0.45 0.45]);
% axis square;
% set(gca,'XTick',[-0.2 -0.1 0 0.1 0.2]);
% set(gca,'YTick',[-0.2 -0.1 0 0.1 0.2]);
% colorbar;
figure;
hold on;
%mesh(xxxx,yyyy,rho,'FaceColor','interp');
mesh(yyyy,xxxx,(log10(rho)),'FaceColor','interp');
%caxis([1.3 3]);
axis([-15 15 -15 15]);
axis square;
print('-depsc','density.eps');
%set(gca,'XTick',[-0.2 -0.1 0 0.1 0.2]);
%set(gca,'YTick',[-0.2 -0.1 0 0.1 0.2]);
colorbar;
%figure;
%spiral(640,10,640,0)
%title('origin');
%figure;
%mesh(yy*dy,xx*dx,(log10(rho)),'FaceColor','interp');
%axis square;
%axis([-15 15 -15 15]);
%title('just rotate')
