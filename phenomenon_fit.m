% reflection of a bloch state by a barrier; numerics vs analytic
% 2016.01.11
clear all; close all; clc; myfont = 22;

step = 0.001;
nPeriod = 11.5;

Llist = [200, 200, 100];
kilist = [80, 100, 50];
Ulist = [1.5, 2, 12];
Tlist = zeros(1, 3);

plist = zeros(6, nPeriod/step + 1);
plist2 = plist;
envolop = zeros(3, nPeriod/step + 1);

for sw = 1: 3
    L = Llist(sw);   N = 2*L+1;
    ki = kilist(sw);
    U = Ulist(sw);
    delta = 2*sin(2*pi/N)*sin(2*pi*ki/N);
    g = U/N;
    T = 2*pi/delta;
    Tlist(sw) = T;
    dt = 0.001*T;
    Tmax = nPeriod*T;
    tlist = 0:dt:Tmax;
    
    theta = 2*atan(g*T);
    omega = theta /T;
    envolop(sw,:) = abs(1 - exp(-i*omega*tlist)).^2/4;
    
    xlist = -L:L;
    xlist = xlist';
    psi0 = (1/sqrt(N))*exp(i*(2*pi*ki/N)*xlist);
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
    
    psi1 = VV'*psi0;
    for s = 1: length(tlist)
        time = (s-1)*dt;
        psi = VV*(exp(-i*time*dd).*psi1);
        plist(1+ (sw-1)*2,s) = abs(psi'*psi0)^2;
        plist(2+ (sw-1)*2,s) = abs(psi'*psif)^2;
        
        n1 = floor(time/T);
        time2 = time - n1*T;
        amp1 = 1/2;
        amp2 = 0.5*(((1-i*2*g*pi/delta)/(1+i*2*g*pi/delta))^n1)*(1-i*2*g*(time2- pi/delta))/(1+ i*2*g*pi/delta);
        plist2(1+ (sw-1)*2,s) = abs( amp1 + amp2 )^2;
        plist2(2+ (sw-1)*2,s) = abs(-amp1 + amp2 )^2;
    end
end

% h1 = figure;
% plot(tlist/T, plist, tlist/T, plist2, ':')
% set(gca,'fontsize',myfont)
% xlim([0 nPeriod])
% ylim([0 1])
% xlabel('t/T','fontsize',myfont);
% ylabel('p','fontsize',myfont);
% str = strcat ('U=', num2str(U),', N=',num2str(N),', ki=',num2str(ki));
% title(str,'fontsize',myfont)

h1= figure;
ha = tight_subplot(1,3,[.05 .05],[.7 .01],[.1 .01])
aa = 0.95;
bb = 0.9;
myfont = 12;

axes(ha(1))
plot((0:step:nPeriod)*Tlist(1), plist(1, :), (0:step:nPeriod)*Tlist(1), plist(2, :),':','linewidth',1.5)
xlabel('$t$','fontsize',myfont, 'Interpreter','Latex')
ylabel('$ P_i \;\& \; P_r $','fontsize',myfont,'Interpreter','Latex')
set(gca,'fontsize',myfont)
xlim([0 nPeriod*Tlist(1)])
ylim([0 1])
XL=xlim; YL=ylim;
str=strcat('(a)');
% text(0.5*XL(1)+0.5*XL(2), 0.5*YL(2)+0.5*YL(1),str,'fontsize',myfont)
text(aa*XL(1)+(1-aa)*XL(2), bb*YL(2)+(1-bb)*YL(1),str,'fontsize',myfont)

% h2=figure;
axes(ha(2))
% set(gca,'position',gcaposition);
plot((0:step:nPeriod)*Tlist(2), plist(3, :), (0:step:nPeriod)*Tlist(2), plist(4, :),':','linewidth',1.5)
xlabel('$t$','fontsize',myfont, 'Interpreter','Latex')
% ylabel('$ |\eta|^2 $','fontsize',myfont,'Interpreter','Latex')
set(gca,'fontsize',myfont)
xlim([0 nPeriod*Tlist(2)])
ylim([0 1])
XL=xlim; YL=ylim;
str=strcat('(b)');
text(aa*XL(1)+(1-aa)*XL(2), bb*YL(2)+(1-bb)*YL(1),str,'fontsize',myfont)

% h3=figure;
% set(gca,'position',gcaposition);
axes(ha(3))
plot((0:step:nPeriod)*Tlist(3), plist(5, :), (0:step:nPeriod)*Tlist(3), plist(6, :),':','linewidth',1.5)
xlabel('$t$','fontsize',myfont, 'Interpreter','Latex')
% ylabel('$ |\eta|^2 $','fontsize',myfont,'Interpreter','Latex')
set(gca,'fontsize',myfont)
xlim([0 nPeriod*Tlist(3)])
ylim([0 1])
XL=xlim; YL=ylim;
str=strcat('(c)');
text(aa*XL(1)+(1-aa)*XL(2), bb*YL(2)+(1-bb)*YL(1),str,'fontsize',myfont)

str = strcat('pipf1.eps');
print(h1,'-depsc',str)

h1= figure;
ha = tight_subplot(1,3,[.05 .05],[.7 .01],[.1 .01])
myfont = 12;

axes(ha(1))
plot((0:step:nPeriod)*Tlist(1), plist(1, :), (0:step:nPeriod)*Tlist(1), plist(2, :),':','linewidth',1.5)
hold on 
plot((0:step:nPeriod)*Tlist(1), plist2(1, :),':', (0:step:nPeriod)*Tlist(1), plist2(2, :),'linewidth',1.5)
% plot((0:step:nPeriod)*Tlist(1),envolop(1,:),'r--','linewidth',1.5)
xlabel('$t$','fontsize',myfont, 'Interpreter','Latex')
ylabel('$ P_i \;\& \; P_r $','fontsize',myfont,'Interpreter','Latex')
set(gca,'fontsize',myfont)
xlim([0 nPeriod*Tlist(1)])
ylim([0 1])
XL=xlim; YL=ylim;
str=strcat('(a)');
% text(0.5*XL(1)+0.5*XL(2), 0.5*YL(2)+0.5*YL(1),str,'fontsize',myfont)
text(aa*XL(1)+(1-aa)*XL(2), bb*YL(2)+(1-bb)*YL(1),str,'fontsize',myfont)

% h2=figure;
axes(ha(2))
% set(gca,'position',gcaposition);
plot((0:step:nPeriod)*Tlist(2), plist(3, :), (0:step:nPeriod)*Tlist(2), plist(4, :),':','linewidth',1.5)
hold on 
plot((0:step:nPeriod)*Tlist(2), plist2(3, :),':', (0:step:nPeriod)*Tlist(2), plist2(4, :),'linewidth',1.5)
% plot((0:step:nPeriod)*Tlist(2),envolop(2,:),'r--','linewidth',1.5)
xlabel('$t$','fontsize',myfont, 'Interpreter','Latex')
% ylabel('$ |\eta|^2 $','fontsize',myfont,'Interpreter','Latex')
set(gca,'fontsize',myfont)
xlim([0 nPeriod*Tlist(2)])
ylim([0 1])
XL=xlim; YL=ylim;
str=strcat('(b)');
text(aa*XL(1)+(1-aa)*XL(2), bb*YL(2)+(1-bb)*YL(1),str,'fontsize',myfont)

% h3=figure;
% set(gca,'position',gcaposition);
axes(ha(3))
plot((0:step:nPeriod)*Tlist(3), plist(5, :), (0:step:nPeriod)*Tlist(3), plist(6, :),':','linewidth',1.5)
hold on 
plot((0:step:nPeriod)*Tlist(3), plist2(5, :),':', (0:step:nPeriod)*Tlist(3), plist2(6, :),'linewidth',1.5)
% plot((0:step:nPeriod)*Tlist(3),envolop(3,:),'r--','linewidth',1.5)
xlabel('$t$','fontsize',myfont, 'Interpreter','Latex')
% ylabel('$ |\eta|^2 $','fontsize',myfont,'Interpreter','Latex')
set(gca,'fontsize',myfont)
xlim([0 nPeriod*Tlist(3)])
ylim([0 1])
XL=xlim; YL=ylim;
str=strcat('(c)');
text(aa*XL(1)+(1-aa)*XL(2), bb*YL(2)+(1-bb)*YL(1),str,'fontsize',myfont)

str = strcat('pipf2.eps');
print(h1,'-depsc',str)