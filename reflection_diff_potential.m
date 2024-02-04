% reflection of a bloch state by a barrier; numerics vs analytic
% 2018.05.07
clear all; close all; clc; myfont = 22;

L = 200;   N = 2*L+1;
ki = 100;
U = 2; lambda = 3;
delta = 2*sin(2*pi/N)*sin(2*pi*ki/N);
g = U/N;
T = 2*pi/delta;
dt = 0.0002*T;
nPeriod = 9.5;
Tmax = nPeriod*T;
tlist = 0:dt:Tmax;
plist = zeros(2, length(tlist));

xlist = -L:L;
xlist = xlist';
psii = (1/sqrt(N))*exp(i*(2*pi*ki/N)*xlist);
psir = (1/sqrt(N))*exp(i*(-2*pi*ki/N)*xlist);
% hamiltonian
H = zeros(N, N);
for s= 1:(N-1)
    H(s,s+1) = -1;     H(s+1,s) = -1;
end
H(1,N) = -1;  H(N,1) = -1;
for s = -L:L
    H(L+1+s, L+1+s) = U* exp(-abs(s)/lambda);
end
[VV,DD] = eig(H);
dd = diag(DD);

psi1 = VV'*psii;
for s = 1: length(tlist)
    time = (s-1)*dt;
    psi = VV*(exp(-i*time*dd).*psi1);
    plist(1,s) = abs(psi'*psii)^2;
    plist(2,s) = abs(psi'*psir)^2;
end

h1 = figure;
plot(tlist/T, plist)
set(gca,'fontsize',myfont)
xlim([0 nPeriod])
ylim([0 1])
xlabel('t/T','fontsize',myfont);
ylabel('p','fontsize',myfont);
str = strcat ('U=', num2str(U),', N=',num2str(N),', ki=',num2str(ki));
title(str,'fontsize',myfont)
