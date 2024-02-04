% spectrum of a toy hamiltonian 
clear all; close all; clc; 

delta = 1;
coupling = 1;

N = 50;

H = zeros(N, N);
for s = 1: N 
    H(s,s) = s*delta;
end

H = H + coupling * ones(N, N);
[VV,DD] = eig(H);
dd = diag(DD);

h1 = figure;
plot(1:N, dd ,'*')

dd(2:end) - dd(1:(end-1))

h2 = figure;
ind =20;
den = abs(VV(:,ind)).^2;
plot(1:N, den, '*' )