% to study the fourier components of S(t)
% 2016.may.06 
clear all; close all; clc; 

delta = 1;
g = 0.2;
L = 60;
N = L + L +1;
T = 2*pi/delta;
theta =2* atan(g*T);
omega = theta / T
steps = 2000;
dt = T/steps;
tlist = dt*(1:steps*2);
Slist = zeros(1, length(tlist));
Slist2 = zeros(1, length(tlist));
amplist = zeros(1, length(tlist));

H = zeros(N, N);
for s = 1: N
    H(s,s) = (s-L-1)*delta;
end

H = H + 2*g * ones(N, N);
[VV,DD] = eig(H);
dd = diag(DD);

amp = zeros(N , 1);
for s1 = 1: N
    amp (s1) = sum(VV(:,s1)) * VV(L+1, s1);
end
max(abs(amp))
amp(end)

sss = (exp(delta/2/g) +1 )/(exp(delta/2/g) -1);
res = delta/(2*g)/(sss/(sss-1) - sss/ (sss+1))

v0 = zeros(N, 1);
v0 (L+1) = 1;
v1 = VV'*v0;
for s = 1: length(tlist)
    v = VV* (exp(-i*dt*(s)*dd).*v1);
    
    Slist(s) = sum(v);
    amplist(s) = v(L+1);    
    Slist2(s) = sum(amp (1:end-1).*exp(-i*dt*s*dd(1:end-1)));
end

h1 = figure;
plot(dd, amp ,'* ')

h2 = figure;
plot3(tlist, real(Slist), imag(Slist))
hold on 
plot3(tlist, real(Slist2), imag(Slist2), 'r')
