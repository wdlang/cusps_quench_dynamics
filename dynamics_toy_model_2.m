% dynamics of the toy model
clear all; close all; clc; myfont = 22;

delta = 1;
g = 0.25;
T = 2*pi/delta;
dt = 0.001*T;
Tmax = 0.5*T;
tlist = 0:dt:Tmax;
plist = zeros(1, length(tlist));

L1 = 200;
L2 = 200;
L = L1+L2+1;
N = 2*L;

H = zeros(N, N);
for s = 1: L 
    H(s,s) = s*delta;
    H(s+L, s+L) = s*delta;
end

H = H + g * ones(N, N);
[VV,DD] = eig(H);
dd = diag(DD);
dd(2:end) - dd(1:(end-1))
% h0 = figure;
% plot(1:N, dd ,'*')
% ylim([0 delta*N/2])
% 
% h2 = figure;
% ind =23;
% den = abs(VV(:,ind)).^2;
% plot(1:N, den, '*' )

v0 = zeros(N, 1);
v0 (L1+1) = 1;

v1 = VV'*v0;

amplist = zeros(1, length(tlist));
plist2 = zeros(1, length(tlist));
plist3 = zeros(1, length(tlist));
interlist = zeros(1, length(tlist));
for s = 1: length(tlist)
    time = dt*(s-1);
    v = VV* (exp(-i*time*dd).*v1);
       
    n1 = floor(time/T);
    time2 = time - n1*T;
    amp1 = 1/2;
    amp2 = 0.5*(((1-i*2*g*pi/delta)/(1+i*2*g*pi/delta))^n1)*(1-i*2*g*(time2- pi/delta))/(1+ i*2*g*pi/delta);
    plist3(s) = abs(amp1 + amp2 )^2;
    plist(s) = abs(v(L1+1))^2;
end

h1 = figure;
plot(tlist/T, plist, tlist/T, plist3,':')
ylim([0 1])
xlabel('t/T','fontsize', myfont)
ylabel('p','fontsize', myfont)
str = strcat ('g= ', num2str(g));
title(str)

% h11 = figure;
% plot3(tlist/T, real(amplist), imag(amplist))
% % ylim([0 1])
% xlabel('t/T','fontsize', myfont)
% ylabel('real','fontsize', myfont)
% zlabel('imaginary','fontsize',myfont)