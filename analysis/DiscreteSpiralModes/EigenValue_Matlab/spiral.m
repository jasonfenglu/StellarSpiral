clear all
close all
GravConst = 4.3e-6;
m = 2; % [number of spiral]
M1=2e10; %[Msun]
M2=5e9; %[Msun]
Ms=5.4e9; %[Msun]
a1 = 12.0; %[kpc]
a2 =8; %[kpc]
as = 1.75; %[kpc]
dr = 0.0001;
r = dr:dr:10.0;
x1 = sqrt(1+(r./a1).^2);
x2 = sqrt(1+(r./a2).^2);
H1 = 59.0625./x1.^11+26.25./x1.^9+16.875./x1.^7+11.25./x1.^5+6.5625./x1.^3;
H2 = 59.0625./x2.^11+26.25./x2.^9+16.875./x2.^7+11.25./x2.^5+6.5625./x2.^3;
sigma0 = 4.5*M1./(pi.*a1^2*x1.^11)-4.5*M2./(pi*a2^2*x2.^11);
Omega = sqrt(16.0/105.0*(M1*GravConst./a1^3.*H1-M2*GravConst./a2^3.*H2)+Ms*GravConst./as^3./(1+(r/as).^2).^1.5);
dOmega_dr = gradient(Omega,dr);
kappa = sqrt(4.0*Omega.^2.*(1+r./(2.0*Omega).*dOmega_dr));
ToomreQ = 1+0.75*exp(-r.^2);
s = -r./Omega.*dOmega_dr;
curF = 2*m*(pi*GravConst*sigma0)./(kappa.^2.*r)./sqrt(1./s-0.5);
snd = ToomreQ.*pi*GravConst.*sigma0./kappa;
omega = 33.360376-1.4113238i;
nu = (omega-m*Omega)./kappa;
k3sqr = (kappa./snd).^2.*(1./ToomreQ.^2-1+nu.^2+0.25*curF.^2.*ToomreQ.^2);
%k3sqr = (kappa./snd).^2.*(1./ToomreQ.^2-1+nu.^2);

temp = abs((real(omega)-m*Omega)./kappa-0.5);
[mintemp indb] = min(temp);
b_old = r(indb);

for k=1:4000
    k
    z = 0;
    u = zeros(size(r));
    u(1) = 1;
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
    
    omega_new=sqrt(-(duoveru(indb)+0.5/sqrt(k3sqr(indb))*dk3_dr(indb))^2/(kappa(indb)/snd(indb))^2-1/(ToomreQ(indb)^2)+1-0.25*curF(indb)^2*ToomreQ(indb)^2)*kappa(indb)+m*Omega(indb)
    %omega_new=sqrt(-(duoveru(indb)+0.5/sqrt(k3sqr(indb))*dk3_dr(indb))^2/(kappa(indb)/snd(indb))^2-1/(ToomreQ(indb)^2)+1)*kappa(indb)+m*Omega(indb)
    
    omega = omega+(omega_new-omega)*0.000001;
    
    nu = (omega-m*Omega)./kappa;
    k3sqr = (kappa./snd).^2.*(1./ToomreQ.^2-1+nu.^2+0.25*curF.^2.*ToomreQ.^2);
    %k3sqr = (kappa./snd).^2.*(1./ToomreQ.^2-1+nu.^2);
    
    temp = abs((real(omega)-m*Omega)./kappa-0.5);
    [mintemp indb] = min(temp);
    b_ne = r(indb)
end  % end k

figure(1)
plot(r,sigma0*1e-6);

figure(2)
plot(r,Omega);

figure(3)
plot(r,Omega);
hold on
plot(r,Omega-kappa/2,'color','r');

figure(4)
plot(r,curF)

figure(5)
plot(r,snd)

figure(6)
plot(r,real(k3sqr))
axis([0 8 -2 2])

figure(7)
plot(r,real(u))
