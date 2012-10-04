%$Id: density.m 112 2011-07-20 03:18:55Z ccfeng $
% THIS ARCHIVE IS GENERATED AUTOMATICALLY BY LAMAC PROGRAM GENERATOR
% WEB INFORMATION http://140.136.192.14/~yen/lamac/lamac2d.php
% Syntax spiral.m
% Version 1.0
% Please do not modify it.
% If you have any problem, please contact with
% lamac@asiaa.sinica.edu.tw
% Date of Generateion :7/21/2004 Time:14:5:34(Hours:Minutes:Seconds)
function density(nend,res)
nbeg = nend;
nstep = 10;
if (nargin<1)
disp_initial
else
disp_evolution(nend,res)
end
function disp_initial
load ini.dat ;
x    =ini(:,1);
y    =ini(:,2);
rho  =ini(:,3);
ux   =ini(:,4);
uy   =ini(:,5);
p    =ini(:,6);
px   =ini(:,7);
py   =ini(:,8);
rsq  = x.^2+y.^2;
r    = sqrt(rsq);
ur   =( ux.*x+uy.*y)./r;
us   =(-ux.*y+uy.*x)./r;
pr   =( px.*x+py.*y)./r;
ps   =(-px.*y+py.*x)./r;
Or   = us./r;
prdOr= pr./r./(Or.^2);
psdOr= ps./r./(Or.^2);
subplot(3,2,1); plot(r,rho)  ; title('density \rho');
subplot(3,2,2); plot(r, ur)  ; title('radial     velocity u_r');
subplot(3,2,3); plot(r, us)  ; title('azmimuthal velocity u_s');
subplot(3,2,4); plot(r,  p)  ; title('potential');
subplot(3,2,5); plot(r,prdOr); title('radial    force/r\Omega^2');
subplot(3,2,6); plot(r,psdOr); title('azimuthal force/r\Omega^2');
function disp_evolution(nend,res)
nbeg = nend;
nstep = 10;
 buffer = 2;
 load variables;
 xmin = variables(1,1); xmax = variables(1,2); ymin = variables(1,3); ymax = variables(1,4);
 pspd = variables(2,1); sspd = variables(2,2); v0sp = variables(2,3); bfam = variables(2,4);
 dcfl = variables(3,1); ncir = variables(3,2); psrl = variables(3,3); ichk = variables(3,4);
 fstr = variables(4,1); fstp = variables(4,2); unk1 = variables(4,3); unk2 = variables(4,4);
 a1   = (nbeg          -mod(nbeg,1000))/1000;
 a2   = (mod(nbeg,1000)-mod(nbeg, 100))/100;
 a3   = (mod(nbeg, 100)-mod(nbeg,  10))/10;
 a4   =  mod(nbeg,10);
 rprt = ['load M' int2str(a1) int2str(a2) int2str(a3) int2str(a4) '.mat'];
 eval(rprt);
 Nx  = sqrt(length(q)/3);
 Ny  = Nx;
 dx  = (xmax-xmin)/(Nx-2*buffer);
 dy  = (ymax-ymin)/(Ny-2*buffer);
 xl  = xmin-dx:dx:xmax+2*dx;
 yl  = ymin-dy:dy:ymax+2*dy;
  x  = xmin-1.5*dx:dx:xmax+1.5*dx;
  y  = ymin-1.5*dy:dy:ymax+1.5*dy;
 [yy,xx]=meshgrid(x,y);
% [xx,yy]=meshgrid(x,y);
 vxs =[xmin,xmax,ymin,ymax];
% RR  = sqrt(xx.^2+yy.^2);
  n = nend;
  a1   = (n          -mod(n,1000))/1000;
  a2   = (mod(n,1000)-mod(n, 100))/100;
  a3   = (mod(n, 100)-mod(n,  10))/10;
  a4   =  mod(n,10);
  rprt = ['load M' int2str(a1) int2str(a2) int2str(a3) int2str(a4) '.mat'];
  eval(rprt);
  q    = reshape(q,Nx,Ny,3);
  rho  = q(:,:,1);
  ux   = q(:,:,2)./rho;
  uy   = q(:,:,3)./rho;
%  ut = sqrt(ux.*ux+uy.*uy);
%  ur   = ( ux.*xx+uy.*yy)./RR;
%  us   = (-ux.*yy+uy.*xx)./RR;
% dyux = ([ux(:,3) ux(:,3:Ny) ux(:,Ny-2)]-[ux(:,1) ux(:,1:Ny-2) ux(:,Ny)])/(2.0*dy);
% dxuy = ([uy(3,:);uy(3:Nx,:);uy(Nx-2,:)]-[uy(1,:);uy(1:Nx-2,:);uy(Nx,:)])/(2.0*dx);
% vort = dxuy-dyux;
  %%% demonstrate the graphs

  l = res * length(q) *2 * pi;
  rr = [1:1:l];
  for i = 1 :l
  	xxx = floor(res*sin(-1*i*2*pi/l)*length(q)/2/max(x))+length(q)/2+1;
  	yyy = floor(res*cos(-1*i*2*pi/l)*length(q)/2/max(y))+length(q)/2+1;
	rr(i)=0.25*(rho(xxx,yyy)+rho(xxx+1,yyy)+rho(xxx,yyy+1)+rho(xxx+1,yyy+1));

  end 
  the = [0:2*180/(l-1):2*180];
  
  plot(the,rr)
  

 



%%%%%%%%%%%%added by jason to save
   FileName = [int2str(n)];
   PathName = ['~/pics/' FileName];
% saveas(gcf,PathName,'jpg')
% print('-depsc',FileName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

   
   
   
   
    title(['\rho(density):', int2str(n),'Myr']);
    rprt = ['print -zbuffer -djpeg ' int2str(a1) int2str(a2) int2str(a3) int2str(a4)];
%   eval(rprt);

%   subplot(2,2,2);
%   mesh(xx,yy,ur,'FaceColor','interp');
%   colorbar;
%   axis(vxs);
%   title('u_r (radial velocity)');
%   subplot(2,2,3);
%   mesh(xx,yy,us,'FaceColor','interp');
%   colorbar;
%   axis(vxs);
%   title('u_s (azimuthal velocity)');
%   subplot(2,2,4);
%   mesh(xx,yy,vort,'FaceColor','interp');
%   colorbar;
%   title('\omega (vorticity)');
%   axis(vxs);
