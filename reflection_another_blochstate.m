% reflection of a bloch state by a barrier
% 2016.01.09
clear all; close all; clc; myfont = 22;

L = 200;   N = 2*L+1;
ki = 100;  increase = -1;
U = 4;
dt = 0.1;     Tmax = 6000;
plist0 = zeros(1, 1+Tmax);
plist = zeros(1, 1+Tmax);
plist_ana = zeros(1, 1+Tmax);
plist_TR = zeros(1, 1+Tmax);

xlist = -L:L;
xlist = xlist';
psi0 = (1/sqrt(N))*exp(i*(2*pi*ki/N)*xlist);
psi00 = (1/sqrt(N))*exp(-i*(2*pi*ki/N)*xlist);
k1 = ki + increase;
psi_another = (1/sqrt(N))*exp(i*(2*pi*k1/N)*xlist);

H = zeros(N, N);
for s= 1:(N-1)
    H(s,s+1) = -1;     H(s+1,s) = -1;
end
H(1,N) = -1;  H(N,1) = -1;
H(L+1, L+1) = U;

[VV,DD] = eig(H);
dd = diag(DD);

plist(1) =0;
psi1 = VV'*psi0;

g = U/N;
Delta = 4*pi*sin(2*pi*ki/N)/N;
T = 2*pi/Delta;
amp = 0.5*(2*g/increase/Delta)*(1/(1+i*g*T));

for s = 1:Tmax
    psi = VV*(exp(-i*dt*s*dd).*psi1);
    plist0 (s+1) = abs(psi'*psi0)^2;
    plist(s+1) = abs(psi'*psi_another)^2;
    plist_ana(s+1) = abs(amp*(1- exp(-i*increase*Delta *dt*s)))^2;
    plist_TR(s+1) = abs(psi'*psi00)^2;
end

h1 = figure;
plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist_ana ,':')
xlabel('t','fontsize',myfont);
ylabel('p','fontsize',myfont);
str = strcat ('U=', num2str(U),', N=',num2str(N),', ki=',num2str(ki));
title(str,'fontsize',myfont)

h2 = figure;
plot(dt*(0:Tmax), plist0+plist_TR)
xlabel('t','fontsize',myfont);
ylabel('p','fontsize',myfont);
str = strcat ('U=', num2str(U),', N=',num2str(N),', ki=',num2str(ki));
title(str,'fontsize',myfont)