% reflection of a bloch state by a barrier; real space; looking for a
% better quantity suitable for experimental observation
% 2016.02.17
clear all; close all; clc; myfont = 22;

L = 200;   N = 2*L+1;
ki = 100;
U = 1.25;
location = 5;
cutoff = 6; 
dt = 0.1;     Tmax = 20000;
plist = zeros(1, 1+Tmax);
plist2 = zeros(1, 1+Tmax);
plist_ana = zeros(1, 1+Tmax);

xlist = -L:L;
xlist = xlist';
psi0 = (1/sqrt(N))*exp(i*(2*pi*ki/N)*xlist);

H = zeros(N, N);
for s= 1:(N-1)
    H(s,s+1) = -1;     H(s+1,s) = -1;
end
H(1,N) = -1;  H(N,1) = -1;
H(L+1, L+1) = U;

[VV,DD] = eig(H);
dd = diag(DD);

plist(1) =1/N ;
psi1 = VV'*psi0;
g = U/N;
Delta = 4*pi*sin(2*pi*ki/N)/N;
T = 2*pi/Delta;
rotation = (1-i*g*T)/(1+i*g*T);

for s = 1:Tmax
    psi = VV*(exp(-i*dt*s*dd).*psi1);   
    plist(s+1) = abs(psi(L+1 + location ))^2;
    
%     amp0=0;
%     for ss = -cutoff: cutoff
%         psi0 = (1/sqrt(N))*exp(i*(2*pi*(ki+ss)/N)*xlist);
%         amp0 = amp0+ (1/sqrt(N))*exp(i*2*pi*(ki+ss)/N*location)*(psi0'*psi);
%         psi0 = (1/sqrt(N))*exp(i*(2*pi*(-ki+ss)/N)*xlist);
%         amp0 = amp0 + (1/sqrt(N))*exp(i*2*pi*(-ki+ss)/N*location)*(psi0'*psi);
%     end
%     plist2(s+1) = N*abs(amp0)^2;
%     
%     p = floor(s*dt/T);
%     tdiff = s*dt - p*T;
%     amp0 = (1 -i*2*g*tdiff/(1+ i*g*T))*(rotation^p);
%     
%     amp = i*sin(2*pi*ki/N*location) + cos(2*pi*ki/N*location) * amp0 ;
%     for ss = 1:cutoff
%         amp = amp + cos(2*pi*(ki+ss)/N*location)* (2*g/ss/Delta )/(1+i *g*T)*(1 - exp(-i*ss*Delta*tdiff))*(rotation^p);
%         amp = amp + cos(2*pi*(ki-ss)/N*location)* (-2*g/ss/Delta )/(1+i *g*T)*(1 - exp(i*ss*Delta*tdiff))*(rotation^p);
%     end
%     
%     plist_ana(s+1) = abs(amp)^2 ;
end
plist = N*plist;

h1 = figure;
plot(dt*(0:Tmax), plist)
set(gca, 'position', [0.15  0.15  0.8  0.8] )
set(gca, 'fontsize', myfont)
% plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist2 ,':')
% plot(dt*(0:Tmax), plist2, dt*(0:Tmax), plist_ana,'--')
xlabel('$t$','fontsize',myfont,'Interpreter','latex');
ylabel('$|\psi_n|^2$','fontsize',myfont,'Interpreter','latex');
% str = strcat ('U=', num2str(U),', N=',num2str(N),', ki=',num2str(ki),', n=',num2str(location));
% title(str,'fontsize',myfont)
% print(h1,'-djpeg','funding4.jpeg')