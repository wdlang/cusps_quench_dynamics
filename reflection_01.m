% reflection of a bloch state by a barrier
% 2016.01.09
clear all; close all; clc; myfont = 22;

L = 200;   N = 2*L+1;
ki = 100;
U = 0.2;
dt = 0.1;     Tmax = 20000;
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

psik = fftshift(fft(psi0));
denk = abs(psik).^2;
h = figure;
hplot= plot(-L:L, denk, '*');
for s = 1:Tmax
    psi = VV*(exp(-i*dt*s*dd).*psi1);
    plist(s+1) = abs(psi'*psi0)^2;

    psik = fftshift(fft(psi));
    denk = abs(psik).^2;
    set(hplot ,'ydata',denk );
    drawnow                     %updates the display
end

h1 = figure;
plot(dt*(0:Tmax), plist)
xlabel('t','fontsize',myfont);
ylabel('p','fontsize',myfont);
str = strcat ('U=', num2str(U),', N=',num2str(N),', ki=',num2str(ki));
title(str,'fontsize',myfont)