% reflection of a bloch state by a barrier; numerics vs analytic
% 2016.01.11
clear all; close all; clc; myfont = 22;

L = 200;   N = 2*L+1;
ki = 100;
U = 0.8;
dt = 0.1;     Tmax = 50000;
plist = zeros(1, 1+Tmax);

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

plist(1) =1;
psi1 = VV'*psi0;

% psik = fftshift(fft(psi0));
% denk = abs(psik).^2;
% h = figure;
% hplot= plot(-L:L, denk, '*');
for s = 1:Tmax
    psi = VV*(exp(-i*dt*s*dd).*psi1);
    plist(s+1) = abs(psi'*psi0)^2;

%     psik = fftshift(fft(psi));
%     denk = abs(psik).^2;
%     set(hplot ,'ydata',denk );
%     drawnow                     %updates the display
end

N1 = 10;
dim = 2*(2*N1+1);
Heff = zeros(dim, dim);
for s = 1:(2*N1+1)
    Heff(s,s) = -2*cos(2*pi*(ki+s-N1-1)/N);
    Heff(s+2*N1+1, s+2*N1+1) = -2*cos(2*pi*(-ki+s-N1-1)/N);
end
Heff = Heff + (U/N)*ones(dim, dim);
v0 = zeros(dim, 1);
v0 (N1+1)= 1;
[VVV,DDD]= eig(Heff);
ddd = diag(DDD);
v0 = VVV'*v0;
plist3 = zeros(1, 1+ Tmax);
for s = 1: length(plist3)
    v = VVV*(exp(-i*dt*(s-1)*ddd).*v0);
    plist3(s) = abs(v(N1+1))^2;
end

h1 = figure;
plot(dt*(0:Tmax), plist, dt*(0:Tmax), plist3, ':')
xlabel('t','fontsize',myfont);
ylabel('p','fontsize',myfont);
str = strcat ('U=', num2str(U),', N=',num2str(N),', ki=',num2str(ki), ', N1=',num2str(N1));
title(str,'fontsize',myfont)