% reflection of a bloch state by a barrier; numerics vs analytic
% 2016.01.11
clear all; close all; clc; myfont = 22;

L = 200;   N = 2*L+1;
ki = 100;
U = 2;
delta = 2*sin(2*pi/N)*sin(2*pi*ki/N);
g = U/N;
T = 2*pi/delta;
dt = 0.0002*T;
nPeriod = 29.5;
Tmax = nPeriod*T;
tlist = 0:dt:Tmax;
plist = zeros(2, length(tlist));

xlist = -L:L;
xlist = xlist';
psii = (1/sqrt(N))*exp(i*(2*pi*ki/N)*xlist);
psif = (1/sqrt(N))*exp(i*(-2*pi*ki/N)*xlist);
% hamiltonian
H = zeros(N, N);
for s= 1:(N-1)
    H(s,s+1) = -1;     H(s+1,s) = -1;
end
H(1,N) = -1;  H(N,1) = -1;
H(L+1, L+1) = U;
[VV,DD] = eig(H);
dd = diag(DD);

psi1 = VV'*psii;
for s = 1: length(tlist)
    time = (s-1)*dt;
    psi = VV*(exp(-i*time*dd).*psi1);
    plist(1,s) = abs(psi'*psii)^2;
    plist(2,s) = abs(psi'*psif)^2;
end

h2 = figure;
plot(tlist, plist(1,:), tlist, plist(2,:),':','linewidth',1.5)
xlim([T - 40, T+40 ])
xlabel('$t$','fontsize',myfont, 'Interpreter','Latex')
ylabel('$ P_i \;\& \; P_r $','fontsize',myfont,'Interpreter','Latex')
set(gca,'fontsize',myfont)

print(h2, '-depsc','details.eps')