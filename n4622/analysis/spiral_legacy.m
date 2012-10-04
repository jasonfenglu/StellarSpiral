%$Id: spiral.m 116 2011-07-20 08:49:39Z ccfeng $
% THIS ARCHIVE IS GENERATED AUTOMATICALLY BY LAMAC PROGRAM GENERATOR
% WEB INFORMATION http://140.136.192.14/~yen/lamac/lamac2d.php
% Syntax spiral.m
% Version 1.0
% Please do not modify it.
% If you have any problem, please contact with
% lamac@asiaa.sinica.edu.tw
% Date of Generateion :7/21/2004 Time:14:5:34(Hours:Minutes:Seconds)
function spiral(nbeg,nstep,nend,ICH)
if (nargin<1)
disp_initial
else
disp_evolution(nbeg,nstep,nend,ICH)
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
function disp_evolution(nbeg,nstep,nend,ICH)
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
 rprt = ['load mats/M' int2str(a1) int2str(a2) int2str(a3) int2str(a4) '.mat'];
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
 for n = nbeg:nstep:nend
  a1   = (n          -mod(n,1000))/1000;
  a2   = (mod(n,1000)-mod(n, 100))/100;
  a3   = (mod(n, 100)-mod(n,  10))/10;
  a4   =  mod(n,10);
  rprt = ['load mats/M' int2str(a1) int2str(a2) int2str(a3) int2str(a4) '.mat'];
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
  if (ICH == 0)
%   subplot(2,2,1);
   mesh(xx,yy,(rho),'FaceColor','interp');
   axis([xmin xmax ymin ymax]); 
%   mesh(xx,yy,rho*153,'FaceColor','interp');
%   contour(xx,yy,ut,20) 
    caxis([8.0  13.0 ]);
%   axis([-5 5 -5 5])
   colorbar;
%   axis(vxs);
%   axis([-0.03 0.03 -0.03 0.03]); 
   axis square;
   tt = (n*790.13)/800;
   title(['\rho(density):', int2str(n)]);


% circles
   hold on
   r = 6.87;
   [sx,sy,sz] = cylinder(r);
   sz = sz /10 ;
%  surf(sx,sy,sz);
   
% circles
   r = 5.76;
   [sx,sy,sz] = cylinder(r);
   sz = sz /10 ;
%  surf(sx,sy,sz);


% circles
   r = 6;
   [sx,sy,sz] = cylinder(r);
   sz = sz /10 ;
   surf(sx,sy,sz);

   hold off
   note = fopen('note','rt');
   comment_position = 0;

%% add comment notes:
   while 1
       comment = fgetl(note);
       if ~ischar(comment), break, end
       text(-30,30-comment_position*2.7,0.05,comment);
       comment_position = comment_position +1;
   end
   fclose(note);
%%%%%%%%%%%%change display range
%  axis([-5 5 -5 5]);



%%%%%%%%%%%%added by jason to save
   FileName = [int2str(n)];
   PathName = ['~/pics/' FileName];
%  saveas(gcf,PathName,'png')
%  print('-depsc',FileName);
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
  elseif (ICH == 1)
   subplot(2,2,1);
   mesh(xx,yy,log10(rho),'FaceColor','interp');
   colorbar;
   axis(vxs);
   title(['\rho(density):', int2str(n)]);
   subplot(2,2,2);
   mesh(xx,yy,ur,'FaceColor','interp');
   colorbar;
   axis(vxs);
   title('u_r (radial velocity)');
   subplot(2,2,3);
   mesh(xx,yy,us,'FaceColor','interp');
   colorbar;
   axis(vxs);
   title('u_s (azimuthal velocity)');
   subplot(2,2,4);
   axis off;
   strcfl   = ['cfl number     =',num2str(dcfl)];
   strSndSpd= ['c(sound speed) =',num2str(sspd)];
   strFrame = ['frame number   =',int2str(n)];
   strGridN = ['Grid Size      =',int2str(Nx-2*buffer)];
   text(-0.1,1.1, strFrame);
   text(-0.1,0.9,   strcfl);
   text(-0.1,0.8,strSndSpd);
  elseif (ICH == 2)
   subplot(2,2,1);
   mesh(xx,yy,log10(rho),'FaceColor','interp');
   colorbar;
   axis(vxs);
   title(['\rho(density):', int2str(n)]);
   subplot(2,2,2);
   mesh(xx,yy,ux,'FaceColor','interp');
   colorbar;
   axis(vxs);
   title('u_x velocity');
   subplot(2,2,3);
   mesh(xx,yy,uy,'FaceColor','interp');
   colorbar;
   axis(vxs);
   title('u_y velocity');
  elseif (ICH == 3)
   subplot(2,2,1);
   contour(xx,yy,rho,20);
   colorbar;
   title(['\rho(density):', int2str(n)]);
   subplot(2,2,2);
   contour(xx,yy,ur,20);
   colorbar;
   title('u_r(radial velocity)');
   subplot(2,2,3);
   contour(xx,yy,us,20);
   colorbar;
   title('u_s(azimuthal velocity)');
  elseif (ICH == 4)
  figure;
   quiver(ux,uy);
  end
  drawnow;
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
