% dynamics of the toy model
clear all; close all; clc; myfont = 22;

delta = 1;
g = 0.125;
T = 2*pi/delta;
steps = 1000;
dt = T/steps;
tlist = dt*(1:steps*30);
plist = zeros(1, length(tlist));
Slist = zeros(1, length(tlist));
Slist2 = zeros(1, length(tlist));

L1 = 10;
L2 = L1;
N = L1 + L2 +1;

H = zeros(N, N);
for s = 1: N
    H(s,s) = (s-L1-1)*delta;
end

H = H + 2*g * ones(N, N);
[VV,DD] = eig(H);
dd = diag(DD);
dd(2:end) - dd(1:(end-1))

v0 = zeros(N, 1);
v0 (L1+1) = 1;

v1 = VV'*v0;
amplist = zeros(1, length(tlist));
amplist2= amplist;
plist2 = amplist;
for s = 1: length(tlist)
    v = VV* (exp(-i*dt*(s)*dd).*v1);
    
    Slist(s) = sum(v);
    amplist(s) = v(L1+1);    
    
    plist2 (s) = abs(1 - i*2*g*dt*(s-1)/(1+i*g*T))^2;
    plist(s) = abs(v(L1+1))^2;
    
    interval = floor((s-1)/steps);
    Slist2(s) = (1/(1+i*g*T))*((1-i*g*T)/(1+i*g*T))^interval;
    
    time = s*dt - interval * T;
    amplist2(s) = ( 1- i*2*g *(time- T/2))/(1+ i*g* T)*((1-i*g*T)/(1+i*g*T))^interval;
end

h1 = figure;
plot(tlist/(2*pi/delta), plist, tlist/(2*pi/delta), plist2,':')
ylim([0 1])
xlabel('t/T','fontsize',myfont);
ylabel('p','fontsize',myfont);
str = strcat ('g=', num2str(g),', N=',num2str(N));
title(str,'fontsize',myfont)

h2 = figure;
plot3(tlist/(2*pi/delta), real(amplist2), imag(amplist2),':','linewidth',1.5,'color','r')
hold on
plot3(tlist/(2*pi/delta), real(amplist), imag(amplist),'linewidth',1.5)
% xlim([0 3])
YL=ylim; ZL=zlim;
xlabel('$t/T $','fontsize', myfont,'Interpreter','latex')
ylabel('$ Re (\psi_0) $','fontsize', myfont,'Interpreter','latex')
zlabel('$ Im (\psi_0 ) $','fontsize', myfont,'Interpreter','latex')
% title('trajectory of \psi_0')
grid on
set(gca,'fontsize', myfont)
set(gca,'LineWidth',2)
YL=ylim; ZL=zlim;
text(0,0.9*(YL(2) -YL(1)) + YL(1),0.9*(ZL(2) - ZL(1))+ZL(1), '(b)','fontsize', myfont)
print(h2,'-depsc','traj_psi0.eps')

h3 = figure;
% plot(tlist/(2*pi/delta), abs(Slist))

hold on 
plot3(tlist(1:steps)/T, real(Slist2(1:steps)), imag(Slist2(1:steps)), '--','color','r','linewidth',1.5)
plot3(tlist((steps+1):2*steps)/T, real(Slist2((steps+1):2*steps)), imag(Slist2((steps+1):2*steps)), '--','color','r','linewidth',1.5)
plot3(tlist(2*steps+1:3*steps)/T, real(Slist2((2*steps+1):3*steps)), imag(Slist2((2*steps+1):3*steps)), '--','color','r','linewidth',1.5)
plot3(tlist/T, real(Slist), imag(Slist),'linewidth',1.5)
view(3)
% xlim([0 3])
xlabel('$t/T $','fontsize', myfont,'Interpreter','latex')
ylabel('$ Re (S) $','fontsize', myfont,'Interpreter','latex')
zlabel('$ Im (S) $','fontsize', myfont,'Interpreter','latex')
grid on
set(gca,'fontsize', myfont)
set(gca,'LineWidth',2)
YL=ylim; ZL=zlim;
text(0,0.9*(YL(2) -YL(1)) + YL(1),0.9*(ZL(2) - ZL(1))+ZL(1), '(a)','fontsize', myfont)
print(h3,'-depsc','traj_S.eps')