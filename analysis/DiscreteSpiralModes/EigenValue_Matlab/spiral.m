clear all
%close all
GravConst = 4.3e-6;
g = GravConst*1e6;
m = 2; % [number of spiral]
dr = 0.01;
r = dr:dr:20.0;


%% Disk Model
Vod = 275.;
rOmega = 2.;
rs = 5.3;
Sod = 1.34;
V = Vod*(1.+rOmega^2./r.^2).^-0.5.*(1.+r.^4./rs^4).^(-0.25*(Sod-1.));
ss = 1206;
h  = 2.83;
sigg = 6.;
sigma0 = ss*(1.-0.45)*exp(-r/h).*f(r/6./h)+6.;

%%%Disk Analysis
%Toomre Q
Qod = 1.0;
q = 2.6;
rq = 3.5;
ToomreQ = Qod*(1.+q*exp(-r.^2/rq^2));

Omega = V./r;

dOmega_dr = gradient(Omega,dr);
kappa = sqrt(4*Omega.^2.*(1+r./(2*Omega).*dOmega_dr));
snd = ToomreQ.*pi*g.*sigma0./kappa;
axis([0 20 0 120]);
s = -r./Omega.*dOmega_dr;
curF = 2*m*(pi*GravConst*sigma0)./(kappa.^2.*r)./sqrt(1./s-0.5);
%omega here==============
omega = 52.220-0.223i;
%========================

nu = (omega-m*Omega)./kappa;
%k3sqr = (kappa./snd).^2.*(1./ToomreQ.^2-1+nu.^2+0.25*curF.^2.*ToomreQ.^2);
k3sqr = (kappa./snd).^2.*(1./ToomreQ.^2-1+nu.^2);

temp = abs((real(omega)-m*Omega)./kappa-0.5);
[mintemp indb] = min(temp);
b_old = r(indb);

for k=1:1
    k
    z = 0;
    u = zeros(size(r));
    u(1) = 10e-7;
    %%%%%%%%% RK4
    for j=2:length(r)
        k1 = z*dr;
        j1 = (-k3sqr(j-1)*u(j-1))*dr;
        %j1 = -u(i-1)*dr;
        
        k2 = (z+0.5*j1)*dr;
        j2 = -(0.5*(k3sqr(j-1)+k3sqr(j)))*(u(j-1)+0.5*k1)*dr;
        %j2 = -(u(i-1)+0.5*k1)*dr;
        
        k3 = (z+0.5*j2)*dr;
        j3 = -(0.5*(k3sqr(j-1)+k3sqr(j)))*(u(j-1)+0.5*k2)*dr;
        %j3 = -(u(i-1)+0.5*k2)*dr;
        
        k4 = (z+j3)*dr;
        j4 = -(k3sqr(j))*(u(j-1)+k3)*dr;
        %j4 = -(u(i-1)+k3)*dr;
        
        u(j) = u(j-1)+(1/6)*(k1+2*k2+2*k3+k4);
        z  =      z+(1/6)*(j1+2*j2+2*j3+j4);
    end
    %%%%%%%%%
    du_dr = gradient(u,dr);
    dk3_dr = gradient(sqrt(k3sqr),dr);
    duoveru=du_dr./u;
    
    err = duoveru(indb)-(-sqrt(k3sqr(indb))*sqrt(-1)-0.5/sqrt(k3sqr(indb))*dk3_dr(indb))
    
    %omega_new=sqrt(-(duoveru(indb)+0.5/sqrt(k3sqr(indb))*dk3_dr(indb))^2/(kappa(indb)/snd(indb))^2-1/(ToomreQ(indb)^2)+1-0.25*curF(indb)^2*ToomreQ(indb)^2)*kappa(indb)+m*Omega(indb)
    omega_new=sqrt(-(duoveru(indb)+0.5/sqrt(k3sqr(indb))*dk3_dr(indb))^2/(kappa(indb)/snd(indb))^2-1/(ToomreQ(indb)^2)+1)*kappa(indb)+m*Omega(indb)
    
    omega = omega+(omega_new-omega)*0.0001;
    
    nu = (omega-m*Omega)./kappa;
    %k3sqr = (kappa./snd).^2.*(1./ToomreQ.^2-1+nu.^2+0.25*curF.^2.*ToomreQ.^2);
    k3sqr = (kappa./snd).^2.*(1./ToomreQ.^2-1+nu.^2);
    
    temp = abs((real(omega)-m*Omega)./kappa-0.5);
    [mintemp indb] = min(temp);
    b_ne = r(indb)

figure(1);
subplot(1,2,1);
plot(r,k3sqr);
title('k3sqr');
axis([0 20 -1 5]);
subplot(1,2,2);
plot(r,u);
title('u');

end  % end k


%figure(1)
%plot(r,sigma0*1e-6);
%
%figure(2)
%plot(r,Omega);
%
%figure(3)
%plot(r,Omega);
%hold on
%plot(r,Omega-kappa/2,'color','r');
%
%figure(4)
%plot(r,curF)
%
%figure(5)
%plot(r,snd)
%
%figure(6)
%plot(r,real(k3sqr))
%axis([0 8 -2 2])
%
%figure(7)
%plot(r,real(u))

