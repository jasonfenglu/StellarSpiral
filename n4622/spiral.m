%$Id: spiral.m 193 2011-09-02 09:12:23Z ccfeng $
% THIS ARCHIVE IS GENERATED AUTOMATICALLY BY LAMAC PROGRAM GENERATOR
% WEB INFORMATION http://140.136.192.14/~yen/lamac/lamac2d.php
% Syntax spiral.m
% Version 1.0
% Please do not modify it.
% If you have any problem, please contact with
% lamac@asiaa.sinica.edu.tw
% Date of Generateion :7/21/2004 Time:14:5:34(Hours:Minutes:Seconds)
function spiral(nbeg,nstep,nend,ICH)
disp_evolution(nbeg,nstep,nend,ICH)
end

function disp_evolution(nbeg,nstep,nend,ICH)
%% loading variables
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
vxs =[xmin,xmax,ymin,ymax];


density = 10*exp(-(xx.^2+yy.^2)*1/48.8581);

%% start playing
for n = nbeg:nstep:nend
    a1   = (n          -mod(n,1000))/1000;
    a2   = (mod(n,1000)-mod(n, 100))/100;
    a3   = (mod(n, 100)-mod(n,  10))/10;
    a4   =  mod(n,10);
    rprt = ['load mats/M' int2str(a1) int2str(a2) int2str(a3) int2str(a4) '.mat'];
    eval(rprt);
    q    = reshape(q,Nx,Ny,3);
    rho  = transpose(q(:,:,1));
    ux   = transpose(q(:,:,2))./rho;
    uy   = transpose(q(:,:,3))./rho;
    
    if (ICH==0)
	    imagesc(x,y,rho);
	    caxis([0 12]);

    elseif (ICH==1)
	    rho  = rho-density;
	    imagesc(x,y,rho);
	    caxis([-0.5 0.5]);
    end


    axis xy;
    axis square;
    colorbar;
    tt = n*2*pi/pspd*ncir/fstp;
    title(['\rho(density):', int2str(n),'  time:',num2str(tt*1000),' Myr']);
    
    
    axis([xmin xmax ymin ymax]);
    
    
    
    %% circles
    %     hold on
    %     r = 6.87;
    %     [sx,sy,sz] = cylinder(r);
    %     sz = sz /10 ;
    %     %  surf(sx,sy,sz);
    
    
    
    %% add comment notes:
    hold off
    note = fopen('note','rt');
    comment_position = 0;
    while 0
        comment = fgetl(note);
        if ~ischar(comment), break, end
        text(-30,30-comment_position*2.7,0.05,comment);
        comment_position = comment_position +1;
    end
    fclose(note);
    
    
    
    %% save image
    FileName = [int2str(n)];
    PathName = ['./pics/' FileName];
       saveas(gcf,PathName,'png')
    %  print('-depsc',FileName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    drawnow;
    pause(0.5);
end
end

